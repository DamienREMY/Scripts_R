
###

rm(list = ls())

working_folder <- "C:/Users/damien.remy/Desktop/global_analysis_podans"

setwd(working_folder)

library(pheatmap)
library(ggplot2)        
library(gridExtra)
library(ggExtra)
library(reshape)
library(ggpubr)

##### Calcul des longueurs de CDS

file_gff <- as.matrix(read.table("Annotation_Mat+_R.gff",
                                 sep = "\t",
                                 header = F))

tab_cds <- file_gff[which(file_gff[,3] == "CDS"),]

taille_tab_cds <- length(tab_cds[,1])

# Création d'un tableau d'annotation avec seulement les CDS
for(i in 1:taille_tab_cds)

{
tab_cds[i,9] <- grep("Pa", tab_cds[i,9], value = TRUE)
tab_cds[i,9] <- sub("gene_id ","\\1",tab_cds[i,9])
tab_cds[i,9]<- sub("transcript.*;","\\1",tab_cds[i,9])
}


# Création d'un tableau vide pour les valeurs de début et de fin de CDS
gene_length <- array(dim = c(taille_tab_cds,3))

# Ajouts des valeurs de début et de fin des CDS
gene_length[,1] <- as.numeric(tab_cds[,"V4"])
gene_length[,2] <- as.numeric(tab_cds[,"V5"])

#Boucle de calcul des longueurs de CDS en kb
for (i in 1:taille_tab_cds){gene_length[i,3] <- (abs(gene_length[i,2] - gene_length[i,1])/1000)}

df_CDS <- data.frame(Taille_CDS=gene_length[,3],
                     ID_CDS=tab_cds[,"V9"],
                     ID_Chr=tab_cds[,"V1"],
                     Sens=tab_cds[,"V7"])

# Addition des lignes de taille respectives de chaque CDS

longueur_CDS <- as.matrix(array(data = tapply(df_CDS$Taille_CDS,df_CDS$ID_CDS,FUN=sum), 
                                dim = c(length(unique(df_CDS$ID_CDS)),1)))

rownames(longueur_CDS) <- unique(df_CDS$ID_CDS)
colnames(longueur_CDS) <- "longueur_CDS"

##### Matrices des htseqcount

## Hisat

htseq_folder_hisat <- paste(working_folder,"/htseq_count_hisat",sep = "")

htseq_files_hisat <- paste("htseq_count_hisat/",
                           grep("hisat_htseqcount.tab",list.files(htseq_folder_hisat),value = TRUE), sep = "")

htseq_cond_hisat <- sub("htseq_count_hisat/","\\1",htseq_files_hisat)
htseq_cond_hisat <- sub("_hisat_htseqcount.tab","\\1",htseq_cond_hisat)

tmp.table <-read.table(htseq_files_hisat[1])

lng.row <- length(longueur_CDS)
lng.col <- length(htseq_files_hisat)

htseq_matrix <- array(dim = c(lng.row,lng.col))

rownames(htseq_matrix) <- tmp.table[1:lng.row,1]
colnames(htseq_matrix) <- htseq_cond_hisat

for (i in 1:lng.col)
{
  tmp.table <-read.table(htseq_files_hisat[i])
  htseq_matrix[,i] <- tmp.table[1:lng.row,2]
  rm(tmp.table)
}

matrix_unorm <-cbind(htseq_matrix,longueur_CDS)

l.row <- length(matrix_unorm[,1])
l.col <- length(matrix_unorm[1,]) - 1

names_col_matrix <- colnames(matrix_unorm[,1:l.col])


names_col_matrix <- gsub(pattern = "ERR22240", x = names_col_matrix, replacement = "Silar")
names_col_matrix <- gsub(pattern = "SRR3197", x = names_col_matrix, replacement = "Lamacchia")
names_col_matrix <- gsub(pattern = "SRR6960", x = names_col_matrix, replacement = "Benocci")
names_col_matrix <- gsub(pattern = "RDBCH", x = names_col_matrix, replacement = "Debuchy")

##### Normalisation des données de comptage


matrix_RPK <- array(data = NA, dim = c(l.row,l.col), dimnames = NULL)
matrix_PM <- array(data = NA, dim = c(1,l.col), dimnames = NULL)
matrix_TPM <- array(data = NA, dim = c(l.row,l.col), dimnames = NULL)

for (i in 1:l.row) 
{
  for (j in 1:l.col)
  { matrix_RPK[i,j] <- (matrix_unorm[i,j] / (matrix_unorm[i,"longueur_CDS"]))
    }}

colnames(matrix_RPK) <- names_col_matrix
rownames(matrix_RPK) <- rownames(matrix_unorm)


write.table(x = matrix_RPK,
            file = "matrix_reads_per_kb.tab",
            sep = "\t",
            quote = F)

for(i in 1:l.col){ matrix_PM[,i] = (sum(matrix_RPK[,i])/1000000)}

write.table(x = matrix_PM,
            file = "matrix_scale_factor.tab",
            sep = "\t",
            quote = F)


for (i in 1:l.row) 
{
  for (j in 1:l.col)
  { matrix_TPM[i,j] <- (matrix_RPK[i,j] / matrix_PM[,j])
    }}

colnames(matrix_TPM) <- colnames(matrix_unorm[,1:l.col])
rownames(matrix_TPM) <- rownames(matrix_unorm)

m.logTPMcount.HISAT <- log10(matrix_TPM + 0.001)

for (i in 1:length(m.logTPMcount.HISAT[1,]))
{
  for (j in 1:length(m.logTPMcount.HISAT[,1]))
  {
    if(m.logTPMcount.HISAT[j,i] < 0)
      {m.logTPMcount.HISAT[j,i] <- 0
    }}}

kenneth.genes <- as.matrix(read.table(file = "kenneth_genes.txt",
                                      header = F))

m.silarTPM <- m.logTPMcount.HISAT[,1:6]

m.kennethTPM <- m.logTPMcount.HISAT[kenneth.genes[,1],1:6]


colnames(m.silarTPM) <- c("M1","M4","P2","P4","A8","A")
colnames(m.kennethTPM) <- c("M1","M4","P2","P4","A8","A")


order.colcount <- c("A","A8","M1","M4","P2","P4")

M4_Genes <- m.silarTPM[,"M4"]

write.table(M4_Genes,
            file = "Genes_M4.tab",
            sep = "\t",
            quote = F)


color_phm <- c("white", "red")

date <- Sys.Date()

fn_pheatmap <- paste("heatmap_clustering_",
                       date,
                       ".png",
                       sep = "")


fn_pheatmap_k <- paste("heatmap_for_kenneth",
                     date,
                     ".png",
                     sep = "")

pheatmap(m.kennethTPM[,order.colcount],
         color = colorRampPalette(color_phm)(10),
         show_rownames = T,
         display_numbers = T,
         cluster_cols = F,
         filename = fn_pheatmap_k)

pheatmap(m.silarTPM[,order.colcount],
         color = colorRampPalette(color_phm)(50),
         show_rownames = F,
         cluster_cols = F,
         filename = fn_pheatmap)

# Matrice avec log10TPM de Silar m.silarTPM

### Histogramme des écarts-types

tab_ecartypeTPM <- array(data = NA, dim = c(l.row,1))
rownames(tab_ecartypeTPM) <- rownames(m.silarTPM)
colnames(tab_ecartypeTPM) <- "ecart.type"

for (i in 1:l.row)
{tab_ecartypeTPM[i,1] <- sd(m.silarTPM[i,])}


df <- data.frame(ecart.type = tab_ecartypeTPM, CDS = rownames(tab_ecartypeTPM))

x11()
par(mfrow=c(1,2))
ftmp <- boxplot(tab_ecartypeTPM,
                main = "Boîte à moustaches des écarts-types respectifs de l'expression de chaque gène au cours du cycle sexuel",
                ylab = "SD log10(TPM)")
hist(tab_ecartypeTPM,
     main = "Histogramme des écarts-types respectifs de l'expression de chaque gène au cours du cycle sexuel",
     xlab = "SD log10(TPM)",
     ylab = "Effectif")

Qt_ecart.type <- ftmp$stats[5,1]

fg1 <- ggplot(df,
       aes(x = ecart.type)) +
  geom_histogram(bins = 500) +
  labs(x = "Ecart-type log10(TPM)",
       y = "Nombre de gènes",
       caption = "")+
  theme(plot.caption = element_text(hjust = 0.5,size=rel(1.2))) +
  coord_cartesian(ylim = c(0,410))+
  geom_vline(xintercept = Qt_ecart.type, linetype = "dashed", colour = "red", size = 1.5)


fg2 <- ggplot(df,
       aes(x = "", y = ecart.type)) +
  geom_boxplot(color = "gray15",
               fill = "ivory3",
               
               # width = 0.1,
               outlier.shape = 16) +
  labs(x = "",
       y = "Ecart-type log10(TPM)",
       caption = "")+
  theme(plot.caption = element_text(hjust = 0.5,size=rel(1.2))) +
  coord_cartesian(ylim = c(0,2.5))

x11()
grid.arrange(fg1, fg2, widths = c(2,1))

### Filtre des gènes 

count_filtered_normLog <- m.silarTPM[which(df[,1] >= Qt_ecart.type),order.colcount]

### K-means

nrow_ratio_km <- c(1:9,seq(from = 10, to = 50, by = 10))

ratio_km <- data.frame(cluster = nrow_ratio_km,
                       ratio = rep(NA,length(nrow_ratio_km)),
                       withinss = rep(NA,length(nrow_ratio_km)),
                       totss = rep(NA,length(nrow_ratio_km)))

i <- 1

for (k in ratio_km[,1])
{
  print(k)
  km_model <- kmeans(count_filtered_normLog,
                     centers = k,
                     nstart = 10,
                     iter.max = 50,
                     algorithm = "Lloyd")
  ratio_km$ratio[i] <- km_model$tot.withinss / km_model$totss
  ratio_km$withinss[i] <- km_model$tot.withinss
  ratio_km$totss[i] <- km_model$totss
  i <- i + 1
}

plot.A <- ggplot(ratio_km, aes(x = cluster)) +
  geom_line(aes(y = ratio, col = "A"),size = 1.2,show.legend = T) +
  labs(title="Représentation graphique de la diminution du ratio distance totale intracluster sur distance totale intercluster \n en fonction du nombre de clusters",
       x="Nombre de clusters",
       y="Distance totale intracluster/distance totale intercluster") +
  guides(color=guide_legend("Légende"))+
  scale_color_manual(name = "colour",
                     labels = "Ratio",
                     values = ("A"="green")) +
  theme(plot.title = element_text(hjust = 0.5)) 

plot.B<- ggplot(ratio_km, aes(x = cluster)) +
  geom_line(aes(y = withinss, col = "A"),size = 1.2,show.legend = T) +
  geom_line(aes(y = totss, col = "B"),size = 1.2,show.legend = T) +
  labs(title="Représentation graphique de la diminution de l'hétérogénéité intracluster \n et de la constance de l'hétérogénéité intercluster en fonction du nombre de clusters",
       x = "Nombre de clusters",
       y = "Distance totale (UA)") +
  guides(color=guide_legend("Légende"))+
  scale_color_manual(name = "colour",
                     labels = c("Distance totale intracluster",
                                "Distance totale intercluster"),
                     values = c("A"="blue",
                                "B"="red"))+  
  theme(plot.title = element_text(hjust = 0.5,
                                  vjust = 1))

x11()
grid.arrange(plot.A,plot.B, ncol = 1, nrow = 2)

centers <- 15

kmeans_countNF <- kmeans(count_filtered_normLog,
                         centers = centers,
                         nstart = 10,
                         iter.max = 50,
                         algorithm = "Lloyd")

meta_cluster <- data.frame("Genes" = kmeans_countNF$size)

fn_km <- paste("heatmap_clustering_kmeans_cycle_",
                     date,
                     ".png",
                     sep = "")

x11()
pheatmap(kmeans_countNF$centers,
         color = colorRampPalette(color_phm)(50),
         show_rownames = T,
         cluster_cols = F,
         display_numbers = T,
         annotation_row = meta_cluster,
         filename = fn_km)

subTPM_KM <- cbind(count_filtered_normLog,kmeans_countNF$cluster)
order_cluster <- order(subTPM_KM[,length(subTPM_KM[1,])]) 
subTPM_KM <- subTPM_KM[order_cluster,]


for (i in 1:centers)
  
{
  gene_list <- which(subTPM_KM[,length(subTPM_KM[1,])] == i)
  gene_list <- rownames(subTPM_KM[gene_list,])
  write(x = gene_list,
        file = paste("genes_list/liste_gene369_cluster_07052019",
                     centers,"_",
                     i,
                     ".txt",
                     sep = "")
        ,sep = "\n")
}








# write.table(m.logTPMcount.HISAT,
#             "logTPMcountHisat.tab",
#             sep = "\t",
#             row.names = T,
#             col.names = T)






#### STAR ##################

htseq_folder_star <- paste(working_folder,"/htseq_count_star",sep = "")

htseq_files_star <- paste("htseq_count_star/",grep("star_htseqcount.tab",list.files(htseq_folder_star),value = TRUE), sep = "")

htseq_cond_star <- sub("htseq_count_star/","\\1",htseq_files_star)
htseq_cond_star <- sub("_star_htseqcount.tab","\\1",htseq_cond_star)

tmp.table <-read.table(htseq_files_star[1])

lng.row <- length(longueur_CDS)
lng.col <- length(htseq_files_star)

htseq_matrix <- array(dim = c(lng.row,lng.col))

rownames(htseq_matrix) <- tmp.table[1:lng.row,1]
colnames(htseq_matrix) <- htseq_cond_star

for (i in 1:lng.col)
{
  tmp.table <-read.table(htseq_files_star[i])
  htseq_matrix[,i] <- tmp.table[1:lng.row,2]
  rm(tmp.table)
}

matrix_unorm <-cbind(htseq_matrix,longueur_CDS)


l.row <- length(matrix_unorm[,1])
l.col <- length(matrix_unorm[1,]) - 1

names_col_matrix <- colnames(matrix_unorm[,1:l.col])
# 
# names_col_matrix <- gsub(pattern = "ERR22240", x = names_col_matrix, replacement = "Silar")
# names_col_matrix <- gsub(pattern = "SRR3197", x = names_col_matrix, replacement = "Lamacchia")
# names_col_matrix <- gsub(pattern = "SRR6960", x = names_col_matrix, replacement = "Benocci")
# names_col_matrix <- gsub(pattern = "RDBCH", x = names_col_matrix, replacement = "Debuchy")

##### Normalisation des données de comptage


matrix_RPK <- array(data = NA, dim = c(l.row,l.col), dimnames = NULL)
matrix_PM <- array(data = NA, dim = c(1,l.col), dimnames = NULL)
matrix_TPM <- array(data = NA, dim = c(l.row,l.col), dimnames = NULL)

for (i in 1:l.row) 
  
{for (j in 1:l.col){ matrix_RPK[i,j] <- (matrix_unorm[i,j] / (matrix_unorm[i,"longueur_CDS"]))}}

colnames(matrix_RPK) <- names_col_matrix
rownames(matrix_RPK) <- rownames(matrix_unorm)

# 
# write.table(x = matrix_RPK,
#             file = "matrix_reads_per_kb.tab",
#             sep = "\t",
#             quote = F)


for(i in 1:l.col){ matrix_PM[,i] = (sum(matrix_RPK[,i])/1000000)}

for (i in 1:l.row) 
{
  for (j in 1:l.col)
  {
    matrix_TPM[i,j] <- (matrix_RPK[i,j] / matrix_PM[,j])
  }}

colnames(matrix_TPM) <- colnames(matrix_unorm[,1:l.col])
rownames(matrix_TPM) <- rownames(matrix_unorm)

m.logTPMcount.STAR <- log10(matrix_TPM + 0.001)

for (i in 1:length(m.logTPMcount.STAR[1,]))
{
  for (j in 1:length(m.logTPMcount.STAR[,1]))
  {
    if(m.logTPMcount.STAR[j,i] < 0){m.logTPMcount.STAR[j,i] <- 0
    }}}


# write.table(matrix_TPM,
#             "TPMcountStar.tab",
#             sep = "\t",
#             row.names = T,
#             col.names = T)



melted_HISAT <- melt(m.logTPMcount.HISAT)

melted_STAR <- melt(m.logTPMcount.STAR)


melt_merge <- cbind(melted_HISAT,melted_STAR[,"value"])
colnames(melt_merge) <- c("Var1","Var2","valueHisat","valueStar")




x11()
ggplot(melt_merge,mapping = aes(x = valueHisat,y = valueStar)) +
geom_jitter()+
  #facet_wrap(. ~ Var2,scales = "free", labeller = label_parsed)+
  stat_cor(method = "spearman", label.x = 2, label.y = 5)+
  labs(caption ="Corrélation des comptages normalisés obtenus entre les aligneurs HISAT et STAR",
       x ="log10(TPM) par gène avec HISAT",
       y ="log10(TPM) par gène avec STAR") +
  xlim(0,6)+
  ylim(0,6)+
  theme(plot.caption = element_text(hjust = 0.5,
                                  vjust = 0.5))

 