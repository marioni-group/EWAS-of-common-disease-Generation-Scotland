#################
### EpiSmokEr ###
#################

# vignette: http://htmlpreview.github.io/?https://github.com/sailalithabollepalli/EpiSmokEr/blob/master/vignettes/epismoker.html

source("http://bioconductor.org/biocLite.R")
library(devtools)
install_github("sailalithabollepalli/EpiSmokEr")

library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(htmlTable)
library(rmarkdown)

suppressPackageStartupMessages({
  library(EpiSmokEr)  
  library(IlluminaHumanMethylation450kmanifest)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(minfi)
  library(htmlTable)
  library(rmarkdown)
})


#### 20k ####

# read in methylation data
GS_w4 <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/mvals.rds")

## CpGs for EpiSmokEr score 
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/epismoker_cpgs.csv")

GS_w4 <- GS_w4[which(row.names(GS_w4) %in% cpgs$CpG),]

# read in sample file
samps=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")

## prepare samplesheet information
samps$Sample_Name <- NULL 
samps$sex <- ifelse(samps$sex=="M", 1, 2)

# match rownames to colnames 
rownames(samps) <- samps$Sample_Sentrix_ID
samps <- samps[rownames(samps) %in% colnames(GS_w4),]
GS_w4 <- GS_w4[, colnames(GS_w4) %in% rownames(samps)]

ids = samps$Sample_Sentrix_ID
GS_w4 <- GS_w4[,match(ids, colnames(GS_w4))]
table(colnames(GS_w4) == rownames(samps))

# convert to beta values
m2beta <- function (M) 
{
  return((2^M)/(2^M + 1))
}

GS_w4 <- m2beta(GS_w4)

# mean impute missing data
# transpose
GS_w4 <- t(GS_w4)

for(i in 1:ncol(GS_w4)){
  GS_w4[is.na(GS_w4[,i]), i] <- mean(GS_w4[,i], na.rm = TRUE)
}

# transpose back
GS_w4 <- t(GS_w4)

## read in files and clean
GS_w4 <- read.csv("U:/Rob/EWAS_Disease_GS/epismoker_dnam.csv")
GS_w4[,1] <- NULL
colnames(GS_w4) <- gsub("X", "", colnames(GS_w4))

t <- read.csv("U:/Rob/EWAS_Disease_GS/epismoker_rownames.csv")
row.names(GS_w4) <- t$V1

samps <- read.csv("U:/Rob/EWAS_Disease_GS/epismoker_samps.csv")
names(samps)[1] <- "Sample_Name"

## Calculate smoking score
result <- epismoker(dataset=GS_w4, samplesheet = samps,method = "SSc")
write.csv(result, "U:/Rob/EWAS_Disease_GS/epismoker_20k.csv",row.names = F)
