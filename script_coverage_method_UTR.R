rm(list=ls())

library(Rcpp)
library(gridExtra)
library(ggplot2)
library(reshape)

working_tolder <- "C:/Users/damien.remy/Desktop/global_analysis_podans"

covfiles_tolder_p <- paste(working_tolder,"/coverage_prom_files_1000",sep = "")
covfiles_tolder_t <- paste(working_tolder,"/coverage_ter_files_1000",sep = "")

setwd(working_tolder)

### Matrice des valeurs de TPM par CDS

matrix_cds_rpk <- as.matrix(read.table("matrix_reads_per_kb.tab",
                                       sep = "\t"))

matrix_scale_factor <- as.matrix(read.table("matrix_scale_factor.tab",
                                 sep = "\t"))






# Récupération des fichiers coverage
prom_files <- paste("coverage_prom_files_1000/",grep("promoters_counts",list.files(covfiles_tolder_p),value = TRUE), sep = "")
prom_cond <- sub(".*coverage_promoters_counts_","\\1",prom_files)
prom_cond <- sub(".txt","\\1",prom_cond)

ter_files <- paste("coverage_ter_files_1000/",grep("terminaters_counts",list.files(covfiles_tolder_t),value = TRUE), sep = "")
ter_cond <- sub(".*coverage_ter_counts_","\\1",ter_files)
ter_cond <- sub(".txt","\\1",ter_cond)


# Nombre de fichiers coverage
n_files <- length(prom_files)

# Initialisation des paramètres des matrices
tmp_matrix_count_p <- as.matrix(read.table(file = prom_files[1],
                                         header = F,
                                         sep = "\t"))

n_area <- as.numeric(length(which(tmp_matrix_count_p[1,9] == tmp_matrix_count_p[,9])))

n_pas <- as.numeric(tmp_matrix_count_p[1,5]) - as.numeric(tmp_matrix_count_p[1,4])

max_area = as.numeric(tmp_matrix_count_p[n_area,5]) - as.numeric(tmp_matrix_count_p[1,4])

# Noms des gènes
row.genes <- tmp_matrix_count_p[,9]
row.genes <- sub(pattern = 'gene_id ',
                 replacement = "\\1",
                 x = row.genes)

row.ugenes <- unique(row.genes)

row.ugenes_ord <- row.ugenes[order(row.ugenes)]

rownames(matrix_cds_rpk) <- row.ugenes_ord




matrix_cds_rpk <- matrix_cds_rpk[row.ugenes,order(colnames(matrix_cds_rpk))]

# Taille des matrices
col.counts <- length(tmp_matrix_count_p[1,])
row.counts <- length(tmp_matrix_count_p[,1])

row.boxplot <- length(row.ugenes)

# Création de la matrice contenant tous les coverage
matrix_counts_p <- as.matrix(array(data = NA, dim = c(row.counts,n_files)))
matrix_counts_t <- as.matrix(array(data = NA, dim = c(row.counts,n_files)))

colnames(matrix_counts_p) <- prom_cond
rownames(matrix_counts_p) <- row.genes

colnames(matrix_counts_t) <- ter_cond
rownames(matrix_counts_t) <- row.genes

# Création de la matrice de somme des coverage
matrix_sum_p <- as.matrix(array(data = NA, dim = c(row.counts,1)))
matrix_sum_t <- as.matrix(array(data = NA, dim = c(row.counts,1)))
# matrix_PM <- array(data = NA, dim = c(1,n_files), dimnames = NULL)

rownames(matrix_sum_p) <- row.genes
rownames(matrix_sum_t) <- row.genes

### Partionnement de la région 5'

start_area <- seq(from = 1, to = max_area, by = n_pas)
end_area <- seq(from = n_pas, to = max_area, by = n_pas)

### Annotation des sous-régions

name_area <- rep(x = NA, times = length(n_area))

# Boucle de l'ajout de chaque fichier coverage dans la matrice matrix_counts_p
for (i in 1:n_files)
  
{
  print(i)
  
  tmp_matrix_count_p <- as.matrix(read.table(file = prom_files[i],
                                           header = F,
                                           sep = "\t"))
  
  tmp_matrix_count_t <- as.matrix(read.table(file = ter_files[i],
                                             header = F,
                                             sep = "\t"))
  
  
  matrix_counts_p[,i] <- tmp_matrix_count_p[,col.counts]
  matrix_counts_t[,i] <- tmp_matrix_count_t[,col.counts]
  
  for (j in 1:row.boxplot)
  {
    
    if (as.numeric(matrix_cds_rpk[j,i]) !=0)
      
    {
      
      k_beg <- 1 + (j-1)*100
      k_end <- 100 + (j-1)*100
      
      
      matrix_counts_p[k_beg:k_end,i] <- (as.numeric(matrix_counts_p[k_beg:k_end,i])/(10*matrix_scale_factor[,i]))/((as.numeric(matrix_cds_rpk[j,i]/(1000*matrix_scale_factor[,i]))))
      matrix_counts_t[k_beg:k_end,i] <- (as.numeric(matrix_counts_t[k_beg:k_end,i])/(10*matrix_scale_factor[,i]))/((as.numeric(matrix_cds_rpk[j,i]/(1000*matrix_scale_factor[,i]))))
      
    }  
      
    else if (as.numeric(matrix_cds_rpk[j,i]) == 0) 
      
    { 
      
      k_beg <- 1 + (j-1)*100
      k_end <- 100 + (j-1)*100
      
      matrix_counts_p[k_beg:k_end,i] <- NA
      matrix_counts_t[k_beg:k_end,i] <- NA
    }
    
   }
  
  # matrix_PM[,i] = (sum(as.numeric(matrix_counts_p[,i])/1000000))
  # matrix_counts_p[,i] <- as.numeric(matrix_counts_p[,i]) / as.numeric(matrix_PM[,i])
  
  rm(tmp_matrix_count_p) # Suppression de la matrice temporaire dans l'environnement
  rm(tmp_matrix_count_t)
}

write.table(x = matrix_counts_p,
            file = "coverage_5UTR.tab",
            sep = "\t")

write.table(x = matrix_counts_t,
            file = "coverage_3UTR.tab",
            sep = "\t")


for(i in 1:row.counts) 
  {

   matrix_sum_p[i] <- median(as.numeric(na.omit(matrix_counts_p[i,])))
   matrix_sum_t[i] <- median(as.numeric(na.omit(matrix_counts_t[i,])))
   
  }

matrix_cov_p <- as.matrix(array(data = NA, dim = c(row.boxplot,n_area)))
matrix_cov_t <- as.matrix(array(data = NA, dim = c(row.boxplot,n_area)))

matrix_sample_cov_p <- as.matrix(array(data = NA, dim = c(row.boxplot,n_area)))
matrix_sample_cov_t <- as.matrix(array(data = NA, dim = c(row.boxplot,n_area)))


rownames(matrix_cov_p) <- row.ugenes
rownames(matrix_cov_t) <- row.ugenes

rownames(matrix_sample_cov_p) <- row.ugenes
rownames(matrix_sample_cov_t) <- row.ugenes

## 13:41 22/05

for(i in 1:length(start_area)) { name_area[i] <- paste(start_area[i],"-",end_area[i],sep = "") }

colnames(matrix_cov_p) <- name_area
colnames(matrix_cov_t) <- name_area

colnames(matrix_sample_cov_p) <- name_area
colnames(matrix_sample_cov_t) <- name_area 

for( i in 1:n_area)
  
{
  tab_part <- seq(from = i, to = row.counts, by = n_area)
  
  matrix_cov_p[,i] <- matrix_sum_p[tab_part]
  matrix_cov_t[,i] <- matrix_sum_t[tab_part]
  
}

for( i in 1:length(matrix_cov_p))
  
{
    matrix_sample_cov_p[i] <- sample(x = matrix_cov_p, size = 1, replace = F)
    matrix_sample_cov_t[i] <- sample(x = matrix_cov_t, size = 1, replace = F)
  
  
}

### Boxplot ggplot2 des données


melted_matrix_cov_p <- melt(matrix_cov_p)
melted_matrix_cov_t <- melt(matrix_cov_t)

melted_sample_cov_p <- melt(matrix_sample_cov_p)
melted_sample_cov_t <- melt(matrix_sample_cov_t)


#melted_matrix_cov <- melt(log10(matrix_cov+0.1))

#melted_sample_cov <- melt(log10(matrix_sample_cov+0.1))



x11()
ggplot(data = melted_matrix_cov_p, aes(y = value)) +
  geom_boxplot(aes (x = X2),
               color = "gray15",
               fill = "ivory3",
               # width = 0.1,
               outlier.shape = 16) +
  coord_cartesian(ylim = c(0,15))+
  scale_x_discrete(breaks = name_area[c(seq(from = 1, to = 100, by = 5),100)],
                   limits = rev(name_area))+
#  scale_y_sqrt()+
  theme(text = element_text(size=20),
        title = element_text(size = 30),
        axis.text.x = element_text(angle=70, hjust=1))+
  xlab("Région en 5' des CDS")+
  ylab("Signal de transcription")


x11()
ggplot(data = melted_matrix_cov_t, aes(y = value)) +
  geom_boxplot(aes (x = X2),
               color = "gray15",
               fill = "ivory3",
               # width = 0.1,
               outlier.shape = 16) +
  coord_cartesian(ylim = c(0,40))+
  scale_x_discrete(breaks = name_area[c(seq(from = 1, to = 100, by = 5),100)],
                   limits = name_area)+
#  scale_y_sqrt()+
  theme(text = element_text(size=20),
        title = element_text(size = 30),
        axis.text.x = element_text(angle=70, hjust=1))+
  xlab("Région en 3' des CDS")+
  ylab("Signal de transcription")


x11()
ggplot(data = melted_sample_cov_p, aes(y = value)) +
  geom_boxplot(aes (x = X2),
               color = "gray15",
               fill = "ivory3",
               # width = 0.1,
               outlier.shape = 16) +
  coord_cartesian(ylim = c(0,10))+
  scale_x_discrete(limits = rev(name_area))

x11()
ggplot(data = melted_sample_cov_t, aes(y = value)) +
  geom_boxplot(aes (x = X2),
               color = "gray15",
               fill = "ivory3",
               # width = 0.1,
               outlier.shape = 16) +
  coord_cartesian(ylim = c(0,40))+
  scale_x_discrete(limits = rev(name_area))

# x11()
# ggplot(df,
#        aes(x = "",y = ecart.type)) +
#   geom_boxplot(color = "gray15",
#                fill = "ivory3",
#                
#                # width = 0.1,
#                outlier.shape = 16) +
#   labs(x = "Ecart-type",
#        y = "Nombre de lectures par CDS",
#        title = "Boîte à moustaches de la distribution des comptages pour chaque gène \n en fonction du stade du cycle sexuel de Podospora anserina")+
#   theme(plot.title = element_text(hjust = 0.5)) +
#   coord_cartesian(ylim = c(0,2.5))







