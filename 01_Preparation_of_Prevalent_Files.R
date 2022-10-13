##################################################
##### EWAS Incident and Prevalent Disease ########
##################################################

## Set working directory 

setwd("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k")

## Load in DNAm data 
start_time <- Sys.time()
meth=readRDS("mvals.rds")
end_time <- Sys.time()
end_time-start_time

## Load in sample information 
samps=readRDS("GS20k_Targets.rds")


###### PREPARE PHENOTYPE FILES FOR PREVALENT DISEASE ########
## Load in disease data 
prev=read.csv("../../GS_dataset/PCQ/disease_19-03-21.csv")


## Merge in Basename information 
prev=merge(prev,samps[,c("Sample_Name", "Sample_Sentrix_ID")],by.x="ID",by.y="Sample_Name", all.y=T)
names(prev)[ncol(prev)] <- "Basename"


## Prepare Alzheimer's proxy phenotype 
prev$AD_Y <- 0 
prev[prev$ID %in% prev[(prev$alzheimers_M %in% 1 | prev$alzheimers_F %in% 1),"ID"],"AD_Y"] <- 1

## Prepare Parkinson's proxy phenotype 
prev$parkinsons_family_Y <- 0
prev[prev$ID %in% prev[(prev$parkinsons_M %in% 1 | prev$parkinsons_F %in% 1),"ID"],"parkinsons_family_Y"] <- 1


## Read in pain information 
pain=read.csv("../../GS_dataset/PCQ/chronic_painv2.csv")
pain1=read.csv("../../GS_dataset/PCQ/chronic_painv5.csv")

## combine v2 and v5 questionnaires and remove duplicates 
pain = rbind(pain[,c("ID","pain_3_months", "back_pain", "neck_pain")],pain1[,c("ID","pain_3_months", "back_pain","neck_pain")])
pain = pain[-which(duplicated(pain$ID)),]

## clean up variables 
pain[pain$back_pain %in% 2,"back_pain"]<-1
pain[pain$neck_pain %in% 2,"neck_pain"]<-1

## merge with prevalent disease file and create pain variable 
prev=merge(prev,pain,by="ID",all.x=T)
prev$pain_Y <- 0
prev[prev$ID %in% pain[(pain$back_pain %in% 1 | pain$neck_pain %in% 1),"ID"],"pain_Y"] <- 1

## Diabetes QC - keep only type 2, remove known type 1 cases or other 
diab.qc=read.csv("../../../Rob/EWAS_Disease_GS/Diabetes.csv")
diab.qc=diab.qc[-which(diab.qc$tname %in% "Type 2"),]
prev[which(prev$ID %in% diab.qc$id),"diabetes_Y"] <- NA

## Create phenotype files 

pheno_ad <- data.frame(FID = prev$Basename, IID = prev$Basename, phen = prev$AD_Y)
write.table(pheno_ad, file="/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/AD.phen", row.names=F, sep=' ')

pheno_park <- data.frame(FID = prev$Basename, IID = prev$Basename, phen = prev$parkinsons_family_Y)
write.table(pheno_park, file="/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/parkinsons_family.phen", row.names=F, sep=' ')


## Breast cancer - just females
pheno_breast <- data.frame(FID = prev$Basename, IID = prev$Basename, phen = prev$breast_cancer_Y)
tmp=samps[samps$sex %in% "F",]
pheno_breast=pheno_breast[pheno_breast$FID %in% tmp$Sample_Sentrix_ID,]
write.table(pheno_breast, file="/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/breast_cancer.phen", row.names=F, sep=' ')

## Prostate cancer - just males 
pheno_prostate <- data.frame(FID = prev$Basename, IID = prev$Basename, phen = prev$prostate_cancer_Y)
tmp=samps[samps$sex %in% "M",]
pheno_prostate=pheno_prostate[pheno_prostate$FID %in% tmp$Sample_Sentrix_ID,]
write.table(pheno_prostate, file="/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/prostate_cancer.phen", row.names=F, sep=' ')

## Remainder as loop 

dis=c("heart_disease_Y", "COPD_Y", "osteo_arthritis_Y", "asthma_Y", "diabetes_Y", "bowel_cancer_Y", "lung_cancer_Y", "stroke_Y", "rheum_arthritis_Y", "pain_Y")

## Loop stage for prevalent phenotypes 
for(i in 1:length(dis)){ 
# get variable 
var=dis[i]

# get phenotype name - for saving file 
phen.name=gsub("_Y", "", var)

# Make phenotype file:
# 3 columns: IID, FID, Phenotype
pheno <- data.frame(FID = prev$Basename, IID = prev$Basename, phen = prev[,as.character(var)])

## save out file 
write.table(pheno, file=paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/", phen.name, ".phen"), row.names=F, sep=' ')
} 

## CKD - added on 29/08/2022
ckd=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/censoring_oct2020_ckd.csv")
pheno <- data.frame(FID = ckd$Sample_Sentrix_ID, IID = ckd$Sample_Sentrix_ID, phen = ckd[,"CKD"])
write.table(pheno, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/CKD.phen", row.names=F, sep=' ')



