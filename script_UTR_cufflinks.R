rm(list = ls())

working_folder <- "C:/Users/Scorvak/Desktop/UTR_analysis"
setwd(working_folder)

library(reshape2)


## Début initialisation

intersect_folder <- paste(working_folder,"/cufflinks",sep = "")

# Récupération des fichiers bedtools
intersect_files <- paste("cufflinks/",grep("merge_",list.files(intersect_folder),value = TRUE), sep = "")
intersect_cond <- sub(".*cufflinks_","\\1",intersect_files)
intersect_cond <- sub("_wb.gff","\\1",intersect_cond)

n_files <- length(intersect_files)


loc_CDS <- as.matrix(read.delim(file = "start_end_CDS_Podans.tab",
                                header = T,
                                sep = "\t"))

n_CDS <- length(loc_CDS[,1])

id_CDS <- loc_CDS[,"ID.Gene"]

matrix_5UTR <- as.matrix(array(data = NA, dim = c(n_CDS,n_files)))
matrix_3UTR <- as.matrix(array(data = NA, dim = c(n_CDS,n_files)))

rownames(matrix_5UTR) <- id_CDS
rownames(matrix_3UTR) <- id_CDS

colnames(matrix_5UTR) <- intersect_cond
colnames(matrix_3UTR) <- intersect_cond

####





## Fin initialisation

for (i in 1:n_files)
  
{
  
  print(i)
  
  tmp_gff <- as.matrix(read.delim(file =intersect_files[i],
                                  header = F,
                                  sep = "\t"))
  
  tmp_gff[,"V5"] <- sub(pattern = ";function.*",
                        replacement = ";",
                        tmp_gff[,"V5"])
  
  tmp_gff[,"V5"] <- sub(pattern = ".*;,",
                        replacement = "\\1",
                        tmp_gff[,"V5"])
  
  tmp_gff[,"V5"] <- sub(pattern = ".*; FPKM",
                        replacement = "FPKM",
                        tmp_gff[,"V5"])
  
  tmp_gff[,"V5"] <- sub(pattern = ".*gene_id",
                        replacement = "gene_id",
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


write.table(matrix_5UTR,
            file = "5UTR_Cufflinks.tab",
            sep = "\t",
            quote = F)

write.table(matrix_3UTR,
            file = "3UTR_Cufflinks.tab",
            sep = "\t",
            quote = F)

matrix_5UTR <- read.table(file = "5UTR_Cufflinks.tab",
                          header = T,
                          sep = " ")

matrix_3UTR <- read.table(file = "3UTR_Cufflinks.tab",
                          header = T,
                          sep = " ")

m_lng5UTR <- array(data = NA,
                   dim = c(n_CDS,n_files))

m_lng3UTR <- array(data = NA,
                   dim = c(n_CDS,n_files))

for(i in 1:n_files)
  
{
  for(j in 1:n_CDS)
    
  {
    m_lng5UTR[j,i] <- abs( as.numeric(matrix_5UTR[j,i]) - as.numeric(loc_CDS[j,"StartCDS"] ))
    
    m_lng3UTR[j,i] <- abs( as.numeric(matrix_3UTR[j,i]) - as.numeric(loc_CDS[j,"EndCDS"]) )
    
  }
  
}


m_med5UTR <- array(data = NA,
                   dim = c(n_CDS,1))

m_med3UTR <- array(data = NA,
                   dim = c(n_CDS,1))

for(i in 1:n_CDS)
  
{
  id_null5 <- which(m_lng5UTR[i,] > 5)
  id_null3 <- which(m_lng3UTR[i,] > 5)
  
  m_med5UTR[i,1] <- median(m_lng5UTR[i,id_null5])
  m_med3UTR[i,1] <- median(m_lng3UTR[i,id_null3])
  
}


melted_lng_5UTR <- melt(m_lng5UTR)
melted_lng_3UTR <- melt(m_lng3UTR)

melted_5UTR <- melt(na.omit(m_med5UTR))
melted_3UTR <- melt(na.omit(m_med3UTR))

write.table(melted_5UTR,
            file = "melted_5UTR_cufflinks.tab",
            sep = "\t",
            quote = F)

write.table(melted_3UTR,
            file = "melted_3UTR_cufflinks.tab",
            sep = "\t",
            quote = F)






