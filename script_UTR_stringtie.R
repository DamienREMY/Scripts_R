rm(list = ls())

library(gridExtra)
library(ggplot2)
library(ggExtra)
library(reshape)

working_folder <- "C:/Users/damien.remy/Desktop/UTR_analysis"
setwd(working_folder)

stringtie_folder <- paste(working_folder,"/stringtie/21052019_c30",sep = "")

stringtie_files <- paste("stringtie/21052019_c30/",grep(".gff",list.files(stringtie_folder),value = TRUE), sep = "")
stringtie_cond <- sub(".*stringtie_","\\1",stringtie_files)
stringtie_cond <- sub("_wb.gff","\\1",stringtie_cond)

n_files <- length(stringtie_files)


loc_CDS <- as.matrix(read.delim(file = "start_end_CDS_Podans.tab",
                                header = T,
                                sep = "\t"))

n_CDS <- length(loc_CDS[,1])

id_CDS <- loc_CDS[,"ID.Gene"]
#
matrix_5UTR <- as.matrix(array(data = NA, dim = c(n_CDS,n_files)))
matrix_3UTR <- as.matrix(array(data = NA, dim = c(n_CDS,n_files)))

rownames(matrix_5UTR) <- id_CDS
rownames(matrix_3UTR) <- id_CDS

colnames(matrix_5UTR) <- stringtie_cond
colnames(matrix_3UTR) <- stringtie_cond

## Fin initialisation


for (i in 1:n_files)

{

  print(i)

  tmp_gff <- as.matrix(read.delim(file = stringtie_files[i],
                                  header = F,
                                  sep = "\t"))

  tmp_gff[,"V5"] <- sub(pattern = ";function.*",
                        replacement = ";",
                        tmp_gff[,"V5"])

  tmp_gff[,"V5"] <- sub(pattern = ".*;,",
                        replacement = "\\1",
                        tmp_gff[,"V5"])

  tmp_gff[,"V5"] <- sub(pattern = "transcript.*;",
                        replacement = "\\1",
                        tmp_gff[,"V5"])

  {
    for (j in 1:n_CDS)

    {

      if (loc_CDS[j,"StrandCDS"] == "+")

      {

        tmp_i <- grep(pattern = loc_CDS[j,"ID.Gene"],
                      x = tmp_gff[,"V5"])

        if ( length(tmp_i) == 1)

        {

          matrix_5UTR[j,i] <- as.numeric(tmp_gff[tmp_i,"V2"])
          matrix_3UTR[j,i] <- as.numeric(tmp_gff[tmp_i,"V3"])

        }

        else if( length(tmp_i) != 1)

        {

          matrix_5UTR[j,i] <- as.numeric(loc_CDS[j,"StartCDS"])
          matrix_3UTR[j,i] <- as.numeric(loc_CDS[j,"EndCDS"])

        }

      }

      if (loc_CDS[j,"StrandCDS"] == "-")

      {

        tmp_i <- grep(pattern = loc_CDS[j,"ID.Gene"],
                      x = tmp_gff[,"V5"])

        if ( length(tmp_i) == 1)

        {
          matrix_5UTR[j,i] <- as.numeric(tmp_gff[tmp_i,"V3"])
          matrix_3UTR[j,i] <- as.numeric(tmp_gff[tmp_i,"V2"])
        }

        else if  ( length(tmp_i) != 1)

        {
          matrix_5UTR[j,i] <- as.numeric(loc_CDS[j,"StartCDS"])
          matrix_3UTR[j,i] <- as.numeric(loc_CDS[j,"EndCDS"])
        }
      }
    }
  }
}


date <- Sys.Date()

fn_5UTR <- paste("matrix_5UTR_stringtie_21052019_c10",
                 date,
                 ".tab",
                 sep = "")

fn_3UTR <- paste("matrix_3UTR_stringtie_21052019_c10",
                 date,
                 ".tab",
                 sep = "")

write.table(matrix_5UTR,
            file = fn_5UTR,
            sep = "\t")

write.table(matrix_3UTR,
            file = fn_3UTR,
            sep = "\t")

matrix_5UTR <- read.table(file = "matrix_5UTR_stringtie_c03062019-05-09.tab",
                          sep = "\t",
                          header = T)

matrix_3UTR <- read.table(file = "matrix_3UTR_stringtie_c03062019-05-09.tab",
                          sep = "\t",
                          header = T) 


m_lng5UTR <- array(data = NA,
                   dim = c(n_CDS,n_files))

m_lng3UTR <- array(data = NA,
                   dim = c(n_CDS,n_files))

rownames(m_lng5UTR) <- id_CDS
rownames(m_lng3UTR) <- id_CDS

colnames(m_lng5UTR) <- stringtie_cond
colnames(m_lng3UTR) <- stringtie_cond 


for(i in 1:n_files)
  
{
  for(j in 1:n_CDS)
    
  {
    
    m_lng5UTR[j,i] <- abs( as.numeric(matrix_5UTR[j,i]) - as.numeric(loc_CDS[j,"StartCDS"]))
    
    m_lng3UTR[j,i] <- abs( as.numeric(matrix_3UTR[j,i]) - as.numeric(loc_CDS[j,"EndCDS"]))
    
  }
  
}

m_med5UTR <- array(data = NA,
                   dim = c(n_CDS,1))

m_med3UTR <- array(data = NA,
                   dim = c(n_CDS,1))


rownames(m_med5UTR) <- id_CDS
rownames(m_med3UTR) <- id_CDS

m_TSS <- array(data = NA, dim = c(n_CDS,1))
m_TES <- array(data = NA, dim = c(n_CDS,1))

colnames(m_TSS) <- "TSS"
colnames(m_TES) <- "TES"

loc_CDS <- cbind(loc_CDS,m_TSS,m_TES)

for(i in 1:n_CDS)
  
{
  id_null5 <- which(m_lng5UTR[i,] > 5)
  id_null3 <- which(m_lng3UTR[i,] > 5)
  
  m_med5UTR[i,1] <- median(m_lng5UTR[i,id_null5])
  m_med3UTR[i,1] <- median(m_lng3UTR[i,id_null3])
  
  if ( is.na(m_med5UTR[i]) == F)
    
  {
    
    if (loc_CDS[i,"StrandCDS"] == "+")
      
    {
      
      loc_CDS[i,"TSS"] <- round(as.numeric(loc_CDS[i,"StartCDS"]) - as.numeric(m_med5UTR[i]))
      loc_CDS[i,"TES"] <- round(as.numeric(loc_CDS[i,"EndCDS"]) + as.numeric(m_med3UTR[i]))
    }
    
    else if (loc_CDS[i,"StrandCDS"] == "-")
      
    {
      loc_CDS[i,"TSS"] <- round(as.numeric(loc_CDS[i,"StartCDS"]) + as.numeric(m_med5UTR[i]))
      loc_CDS[i,"TES"] <- round(as.numeric(loc_CDS[i,"EndCDS"]) - as.numeric(m_med3UTR[i]))}
    
  }
  
  if (is.na(loc_CDS[i,"TSS"]) == T)
    
  {
    loc_CDS[i,"TSS"] <- as.numeric(loc_CDS[i,"StartCDS"])
    
  }
  
  if (is.na(loc_CDS[i,"TES"]) == T)
    
  {
    loc_CDS[i,"TES"] <- as.numeric(loc_CDS[i,"EndCDS"])
    
  }
  
}
  




write.table(x = loc_CDS,
            file = "TSS_Genes_Podans_c10.tab",
            sep = "\t",
            quote = F)





### Synthèse du fichier .bed des TSS

M4_Genes <- read.table(file = "Genes_M4.tab",
                       header = F,
                       sep = "\t")


ExpG <- paste(M4_Genes[which(M4_Genes[,2] > 2),1],";",sep = "")
UnxpG <- paste(M4_Genes[which(M4_Genes[,2] < 0.5),1],";",sep = "")

Exp_TSS_TES <- array(data = NA, dim = c(length(ExpG),3))
Unxp_TSS_TES <- array(data = NA, dim = c(length(UnxpG),3))

for(i in 1:length(ExpG))

{
  Exp_TSS_TES[i,] <- loc_CDS[grepl(pattern = ExpG[i], x = loc_CDS[,"ID.Gene"], fixed = T),6:8]
  
}
 
for(i in 1:length(UnxpG))
  
{
  Unxp_TSS_TES[i,] <- loc_CDS[grepl(pattern = UnxpG[i], x = loc_CDS[,"ID.Gene"], fixed = T),6:8]
  
}

write.table(Exp_TSS_TES,
            file = "bed_Exp_TSS_TES_21052019.bed",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = F)




write.table(Unxp_TSS_TES,
            file = "bed_Unxp_TSS_TES_21052019.bed",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = F)



melted_lng_5UTR <- melt(m_lng5UTR)
melted_lng_3UTR <- melt(m_lng3UTR)

melted_median_5UTR <- melt(na.omit(m_med5UTR))
melted_median_3UTR <- melt(na.omit(m_med3UTR))

# erreur mean_5stringtie <- mean(melted_median_5UTR[,"value"])
# erreur mean_5stringtie <- mean(m_med5UTR[,"value"])

mean_5stringtie <- median(melted_median_5UTR[,"value"])

fg5UTR_stringtie <- ggplot(melted_median_5UTR,
                           aes(x = value)) +
  geom_histogram(color = "gray15",
               fill = "ivory3",
               binwidth = 50) +
  scale_x_continuous(name = "Taille des régions 5'UTR (pb) par Stringtie",
                     breaks = seq(from = 0, to = 3000, by = 250))+
  labs(y = "Nombre de gènes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = mean_5stringtie, linetype = "dashed", colour = "red", size = 1.5)+
  coord_cartesian(xlim = c(0,3000)
                  ,ylim = c(0,1500))

mean_3stringtie <- median(melted_median_3UTR[,"value"])

fg3UTR_stringtie <- ggplot(melted_median_3UTR,
                           aes(x = value)) +
  geom_histogram(color = "gray15",
               fill = "ivory3",
               binwidth = 50) +
  scale_x_continuous(name = "Taille des régions 3'UTR (pb) par Stringtie",
                     breaks = seq(from = 0, to = 3000, by = 250))+
  labs(y = "Nombre de gènes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = mean_3stringtie, linetype = "dashed", colour = "red", size = 1.5)+
  coord_cartesian(xlim = c(0,3000)
                  ,ylim = c(0,1500))






# 
# x11()
# ggplot(melted_lng_5UTR,
#        aes(x = value))+
#   geom_histogram(binwidth = 50)+
#   facet_wrap(.~ X2, scales = "free", labeller = label_parsed)+
#   theme(axis.title.y =element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.x=element_blank())+
#   coord_cartesian(xlim = c(0,2000),ylim = c(0,500))+
#   labs(x = "RNA-seq",
#        y = "Taille des régions 5'UTR",
#        title = " Histogrammes de distribution des tailles des régions 5'UTR obtenues par Stringtie \n en fonction des RNA-seq analysés")
#   
# x11()
# ggplot(melted_lng_3UTR,
#        aes(x = value))+
#   geom_histogram(binwidth = 50)+
#   facet_wrap(.~ X2, scales = "free", labeller = label_parsed)+
#   theme(axis.title.y =element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.x=element_blank())+
#   coord_cartesian(xlim = c(0,2000),ylim = c(0,500))+
#   labs(x = "RNA-seq",
#        y = "Taille des régions 3'UTR",
#        title = " Histogrammes de distribution des tailles des régions 3'UTR obtenues par Stringtie \n en fonction des RNA-seq analysés")

#### Cufflinks 

melted_5UTR_cufflinks <- read.table(file = "melted_5UTR_cufflinks.tab",
                                    sep = "\t",
                                    header = T)

#mean_5cufflinks <- mean(melted_5UTR_cufflinks[,"value"])
mean_5cufflinks <- median(melted_5UTR_cufflinks[,"value"])

fg5UTR_cufflinks <- ggplot(melted_5UTR_cufflinks,
                           aes(x = value)) +
  geom_histogram(color = "gray15",
                 fill = "ivory3",
                 binwidth = 50) +
  scale_x_continuous(name = "Taille des régions 5'UTR (pb) par Cufflinks",
                     breaks = seq(from = 0, to = 3000, by = 250))+
  labs(y = "Nombre de gènes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = mean_5cufflinks, linetype = "dashed", colour = "red", size = 1.5)+
  coord_cartesian(xlim = c(0,3000)
                  ,ylim = c(0,1500))



melted_3UTR_cufflinks <- read.table(file = "melted_3UTR_cufflinks.tab",
                                    sep = "\t",
                                    header = T)

#mean_5cufflinks <- mean(melted_3UTR_cufflinks[,"value"])

mean_3cufflinks <- median(melted_3UTR_cufflinks[,"value"])



fg3UTR_cufflinks <- ggplot(melted_3UTR_cufflinks,
                           aes(x = value)) +
  geom_histogram(color = "gray15",
                 fill = "ivory3",
                 binwidth = 50) +
  scale_x_continuous(name = "Taille des régions 3'UTR (pb) par Cufflinks",
                     breaks = seq(from = 0, to = 3000, by = 250))+
  labs(y = "Nombre de gènes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = mean_3cufflinks, linetype = "dashed", colour = "red", size = 1.5)+
  coord_cartesian(xlim = c(0,3000)
                  ,ylim = c(0,1500))


x11()
grid.arrange(fg5UTR_stringtie,
             fg5UTR_cufflinks,
             ncol = 1,
             nrow = 2)

x11()
grid.arrange(fg3UTR_stringtie,
             fg3UTR_cufflinks,
             ncol = 1,
             nrow = 2)






