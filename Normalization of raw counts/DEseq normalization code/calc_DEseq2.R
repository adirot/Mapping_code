#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()

#BiocManager::install(c("htmltools","DESeq2"),force = TRUE)
#install.packages("htmltools")
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

library(DESeq2)

#Experiment_name = "Experiment1"
#Experiment_name = "Experiment2"
#comp_strs <- list('CASP','Disrupted')

Experiment_name = "Exponential"
comp_strs <- list('Exponential')


for (comp_str in comp_strs){
  countData <- read.csv(paste(getwd( ),'/raw_data_organized_for_R_',Experiment_name,'/',comp_str,'_raw_reads_for_R.csv',sep=""), header = TRUE, sep = ",")
  head(countData)
  
  
  metaData <- read.csv(paste(getwd( ),'/raw_data_organized_for_R_',Experiment_name,'/',comp_str,'_metadata_for_R.csv',sep=""), header = TRUE, sep = ",")
  metaData
  
  # create DEseq2 object
  dds <- DESeqDataSetFromMatrix(countData=countData, 
                                colData=metaData, 
                                design=~ 1, tidy = TRUE)
  
  
  # perform DEseq2 
  dds <- DESeq(dds)
  
  if (!dir.exists(paste(getwd( ),'/DEseq2_normalized_',Experiment_name,sep=""))){
    dir.create(paste(getwd( ),'/DEseq2_normalized_',Experiment_name,sep=""))
  }
  
  write.csv(counts(dds, normalized=T),paste(getwd( ),'/DEseq2_normalized_',Experiment_name,'/',comp_str,'_DEseq2_norm.csv',sep=""))
}