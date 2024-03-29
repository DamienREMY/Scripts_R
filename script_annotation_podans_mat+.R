
# Script �crit par Damien REMY, contact : damien.remy@u-psud.fr

### Fichier � avoir dans le dossier de script:

# - fichier d'annotation du g�nome de Podospora anserina S mat + 'annotationPodo.gff'
# - liste des fonctions par g�ne 'function_genes_Podo_mat+.txt'

# Que fait ce script ? 

# Ce script modifie et compl�te la derni�re colonne (9) du fichier d'annotation 'annotationPodo.gff'

##### ##### ##### ##### ##### 

rm(list=ls())

htseq_folder <- 'C:/Users/damien.remy/Desktop/Script_annotation_podans'
setwd(htseq_folder)

# R�cup�ration du fichier d'annotation du g�nome de PODANS
GTF_Podo <- as.matrix(read.delim('annotationPodo.gff',header = FALSE,sep = '\t'))

# R�cup�ration de la liste des fonctions par g�ne de PODANS (type sexuel mat + seulement)
Function_Podo <- as.matrix(read.delim('function_genes_Podo_mat+.txt',
                                      header = F,
                                      sep = '\t'))

# Extraction de l'annotation du g�nome mat+

GTF_Podo_P <- GTF_Podo[which(GTF_Podo[,1] != 'Chr1moins'),]

GTF_Podo_P[,4] <- as.numeric(GTF_Podo_P[,4])
GTF_Podo_P[,5] <- as.numeric(GTF_Podo_P[,5])

# Colonne ID g�nes, exons

ID_plus <- GTF_Podo_P[,'V9']
lnp <- length(ID_plus)

ex <- replicate(lnp,1)

# Boucle num�rotation des exons pour mat+

for (i in 1:(lnp-1))
  
{
  if (grepl(fixed = TRUE,
            pattern = paste(ID_plus[i],';',sep = ''),
            x = paste(ID_plus[i+1],';', sep = '')))

  {ex[i+1] <- ex[i] + 1}}

ID_tmp <- paste(ID_plus,';',sep = '')

#D�but de boucle pour compl�ter le fichier d'annotation du type sexuel mat+

for (i in 1:lnp)

{
  if (substr(ID_plus[i],1,2) == 'Pa')
    
  { # Annotation des CDS
    
    ID_plus[i] <- paste('gene_id "',ID_plus[i],'";'
                        ,'transcript_id "',ID_plus[i],'";'
                        ,'exon_number "',ex[i],'";'
                        ,'function "',Function_Podo[which(ID_tmp[i] == paste(Function_Podo[,1],';',sep = '')),2],'";'
                        ,'color "#0000EE";'
                        , sep = '')}
  
  else if (grepl(fixed = TRUE,GTF_Podo_P[i,'V3'],'repeat_region') | grepl(fixed = TRUE,GTF_Podo_P[i,'V3'],'repeat_unit'))
  
  { # Annotation des s�quences r�p�t�es
    
    ID_plus[i] <- paste('id "repeat_sequence_',ID_plus[i],'";'
                        ,'color "#800000";'
                        , sep = '')}
  
  else if (grepl(fixed = TRUE,GTF_Podo_P[i,'V3'],'tRNA'))
    
  { # Annotation des tRNA
    
    ID_plus[i] <- paste('id "tRNA','";'
                        ,'color "#a5a5ff";'
                        , sep = '')}
  
  else if (grepl(fixed = TRUE,GTF_Podo_P[i,'V3'],'rRNA'))
    
  { # Annotation des rRNA
    
    ID_plus[i] <- paste('id "rRNA:','";'
                        ,'color "#cebc96";'
                        , sep = '')}} #Fin de boucle annotation Podo S mat+

# Remplacement de la colonne d'annotation

GTF_Podo_P[,'V9'] <- ID_plus

tab_cds <- GTF_Podo_P[which(GTF_Podo_P[,3] == "CDS"),]

write.table(GTF_Podo_P,
            'Annotation_Mat+_R.gff',
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t',
            quote = F)

write.table(tab_cds,
            'Annotation_Mat+_CDS_R.gff',
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t',
            quote = F)
