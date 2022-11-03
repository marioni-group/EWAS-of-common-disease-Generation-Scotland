## Load requisite libraries 
library(data.table)
library(coloc)
## Read in CpGs for testing 
## Incident disease ##
cpgs=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/incident_full.csv")
cpgs=cpgs[which(cpgs$Present_in_Basic_Model==1),]
# Subset to GoDMC files 
files=list.files("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/mQTLs/GoDMC/", ".")
files=gsub(".fastGWA.rsid", "", files)
cpgs=cpgs[cpgs$Probe %in% files,]
# Set up summary dataframe 
output <- matrix(nrow = nrow(cpgs), ncol = 8)
output <- as.data.frame(output)
names(output) <- c("CpG", "Trait", "nsnps", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")
# Loop through CpGs 
for(i in 1:nrow(cpgs)){
# subset dataframe to just that cpg 
tmp=cpgs[i,]
# extract trait 
trait=tmp$trait
# extract chr and position 
chr=tmp$Chr
pos=tmp$bp 
pos1=tmp$bp+5e5
pos2=tmp$bp-5e5
# read in summary statistics for this dataframe 
tmp1=as.data.frame(fread(paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/mQTLs/GoDMC/", as.character(tmp[1,"Probe"]), ".fastGWA.rsid")))
tmp2=tmp1[which(tmp1$CHR==chr),]
tmp2=tmp2[which(tmp2$POS <= pos1),]
tmp2=tmp2[which(tmp2$POS >= pos2),]
if(nrow(tmp2)==0){ NULL } else {
# read in associated disease state 
dis=as.data.frame(fread(paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/GWAS/Cleaned/", trait, ".txt")))
dis2=dis[which(dis$CHR==chr),]
dis2=dis2[which(dis2$POS <= pos1),]
dis2=dis2[which(dis2$POS >= pos2),]
# create MAF columns 
tmp2$MAF <- ifelse(tmp2$FREQ > 0.5, 1-tmp2$FREQ, tmp2$FREQ)
dis2$MAF <- ifelse(dis2$FREQ > 0.5, 1-dis2$FREQ, dis2$FREQ)
# clean dataframes prior to colocalisation 
## Missing MAF 
if(length(which(is.na(dis2$MAF))) >= 1) { 
  dis2 <- dis2[-which(is.na(dis2$MAF)),]
}
if(length(which(is.na(tmp2$MAF))) >= 1) { 
  tmp2 <- tmp2[-which(is.na(tmp2$MAF)),]
} 
## Duplicated SNPs 
if(length(which(duplicated(dis2$SNP))) >= 1) { 
  dis2 <- dis2[-which(duplicated(dis2$SNP)),]
} 
if(length(which(duplicated(tmp2$rsid))) >= 1) { 
  tmp2 <- tmp2[-which(duplicated(tmp2$rsid)),]
} 
## Missing SNPs
if(length(which(is.na(dis2$SNP))) >= 1) { 
  dis2 <- dis2[-which(is.na(dis2$SNP)),]
} 
if(length(which(is.na(tmp2$rsid))) >= 1) { 
  tmp2 <- tmp2[-which(is.na(tmp2$rsid)),]
} 
# prepare dataframes for coloc 
dataset1 = list(snp = as.character(tmp2$rsid), N = as.numeric(tmp2$N), beta = as.numeric(tmp2$BETA), varbeta = as.numeric(tmp2$SE)^2, MAF = as.numeric(tmp2$MAF), type = "quant")
dataset2 = list(snp = as.character(dis2$SNP), N = as.numeric(dis2$N), beta= as.numeric(dis2$BETA), MAF = as.numeric(dis2$MAF), varbeta = as.numeric(dis2$SE)^2, type = "cc")
# run coloc 
coloc = coloc.abf(dataset2, dataset1)
# store output 
ind = i 
output[ind, 1] <- as.character(tmp[1,"Probe"])
output[ind, 2] <- as.character(trait)
output[ind, 3] <- coloc[[1]][[1]]
output[ind, 4] <- coloc[[1]][[2]]
output[ind, 5] <- coloc[[1]][[3]]
output[ind, 6] <- coloc[[1]][[4]]
output[ind, 7] <- coloc[[1]][[5]]
output[ind, 8] <- coloc[[1]][[6]]
# print to denote completion 
print(ind)
} 
}
# Save out file for incident disease 
output=output[which(!is.na(output$CpG)),]
write.csv(output, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/incident_coloc_godmc.csv",row.names=F)


## Prevalent disease ##
cpgs=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/prevalent_full.csv")
cpgs=cpgs[which(cpgs$Present_in_Basic_Model==1),]
# Subset to GoDMC files 
files=list.files("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/mQTLs/GoDMC/", ".")
files=gsub(".fastGWA.rsid", "", files)
cpgs=cpgs[cpgs$Probe %in% files,]
# Set up summary dataframe 
output <- matrix(nrow = nrow(cpgs), ncol = 8)
output <- as.data.frame(output)
names(output) <- c("CpG", "Trait", "nsnps", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")
# Loop through CpGs 
for(i in 1:nrow(cpgs)){
  # subset dataframe to just that cpg 
  tmp=cpgs[i,]
  # extract trait 
  trait=tmp$trait
  # extract chr and position 
  chr=tmp$Chr
  pos=tmp$bp 
  pos1=tmp$bp+5e5
  pos2=tmp$bp-5e5
  # read in summary statistics for this dataframe 
  tmp1=as.data.frame(fread(paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/mQTLs/GoDMC/", as.character(tmp[1,"Probe"]), ".fastGWA.rsid")))
  tmp2=tmp1[which(tmp1$CHR==chr),]
  tmp2=tmp2[which(tmp2$POS <= pos1),]
  tmp2=tmp2[which(tmp2$POS >= pos2),]
  if(nrow(tmp2)==0){ NULL } else { 
  # read in associated disease state 
  dis=as.data.frame(fread(paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/GWAS/Cleaned/", trait, ".txt")))
  dis2=dis[which(dis$CHR==chr),]
  dis2=dis2[which(dis2$POS <= pos1),]
  dis2=dis2[which(dis2$POS >= pos2),]
  # create MAF columns 
  tmp2$MAF <- ifelse(tmp2$FREQ > 0.5, 1-tmp2$FREQ, tmp2$FREQ)
  dis2$MAF <- ifelse(dis2$FREQ > 0.5, 1-dis2$FREQ, dis2$FREQ)
  # clean dataframes prior to colocalisation 
  ## Missing MAF 
  if(length(which(is.na(dis2$MAF))) >= 1) { 
    dis2 <- dis2[-which(is.na(dis2$MAF)),]
  }
  if(length(which(is.na(tmp2$MAF))) >= 1) { 
    tmp2 <- tmp2[-which(is.na(tmp2$MAF)),]
  } 
  ## Duplicated SNPs 
  if(length(which(duplicated(dis2$SNP))) >= 1) { 
    dis2 <- dis2[-which(duplicated(dis2$SNP)),]
  } 
  if(length(which(duplicated(tmp2$rsid))) >= 1) { 
    tmp2 <- tmp2[-which(duplicated(tmp2$rsid)),]
  } 
  ## Missing SNPs
  if(length(which(is.na(dis2$SNP))) >= 1) { 
    dis2 <- dis2[-which(is.na(dis2$SNP)),]
  } 
  if(length(which(is.na(tmp2$rsid))) >= 1) { 
    tmp2 <- tmp2[-which(is.na(tmp2$rsid)),]
  } 
  # prepare dataframes for coloc 
  dataset1 = list(snp = as.character(tmp2$rsid), N = as.numeric(tmp2$N), beta = as.numeric(tmp2$BETA), varbeta = as.numeric(tmp2$SE)^2, MAF = as.numeric(tmp2$MAF), type = "quant")
  dataset2 = list(snp = as.character(dis2$SNP), N = as.numeric(dis2$N), beta= as.numeric(dis2$BETA), MAF = as.numeric(dis2$MAF), varbeta = as.numeric(dis2$SE)^2, type = "cc")
  # run coloc 
  coloc = coloc.abf(dataset2, dataset1)
  # store output 
  ind = i
  output[ind, 1] <- as.character(as.character(tmp[1,"Probe"]))
  output[ind, 2] <- as.character(trait)
  output[ind, 3] <- coloc[[1]][[1]]
  output[ind, 4] <- coloc[[1]][[2]]
  output[ind, 5] <- coloc[[1]][[3]]
  output[ind, 6] <- coloc[[1]][[4]]
  output[ind, 7] <- coloc[[1]][[5]]
  output[ind, 8] <- coloc[[1]][[6]]
  # print to denote completion 
  print(ind)
} 
}
output=output[which(!is.na(output$CpG)),]
write.csv(output, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/prevalent_coloc_godmc.csv",row.names=F)
