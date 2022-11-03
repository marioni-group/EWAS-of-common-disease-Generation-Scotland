## Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k")

## Load in sample information 
samps=readRDS("GS20k_Targets.rds")

## Subset to covariates of interest 
samps1 = samps[,c("Sample_Name", "Sample_Sentrix_ID", "age", "sex", "Set", "Batch")]

## Read in additional covariates for fully-adjusted model 
# pack years 
smk = read.csv("../../GS_dataset/updated_smoking_jan_2019/pack_years.csv")
# bmi 
bmi = read.csv("../../GS_dataset/clinical/body.csv")
names(bmi)[1] <- "Sample_Name"
bmi <- bmi[,c("Sample_Name", "bmi")]
# cell counts 
w1 = readRDS("../stradl-samples-5087.rds")
w1 <- w1[,c("Sample_Name", "Bcell", "CD4T", "CD8T", "Gran", "NK", "Mono")]

w3 = read.csv("../wave3-final/samplesheet.final.csv")
w3 <- w3[,c("Sample_Name", "Bcell", "CD4T", "CD8T", "Gran", "NK", "Mono")]

w4 = read.table("../wave4/wave4_cellcomp.tsv",header = T)
w4 = merge(w4, samps[,c("Sample_Sentrix_ID", "Sample_Name")], by.x = "ID", by.y = "Sample_Sentrix_ID")
w4$ID <- NULL 
w4 <- w4[,c("Sample_Name","Bcell", "CD4T", "CD8T", "Gran", "NK", "Mono")]

cells = rbind(w1, rbind(w3, w4))

## simd 
simd=read.csv("../../GS_dataset/clinical/SIMD.csv")
names(simd)[1] <- "Sample_Name"
simd<-simd[,c("Sample_Name", "rank")]

## EA 
ea=read.csv("../../GS_dataset/PCQ/education.csv")
names(ea)[1] <- "Sample_Name"
ea<-ea[,c("Sample_Name", "years")]

## Alc 
alc = read.csv("../../GS_dataset/PCQ/alcohol.csv")
names(alc)[1] <- "Sample_Name"
alc <- alc[,c("Sample_Name", "units", "usual")]

## epismoker score 
epismk = read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/epismoker_20k.csv")
epismk = merge(epismk, samps[,c("Sample_Sentrix_ID", "Sample_Name", "Set", "Batch", "age", "sex")], by.x = "SampleName", by.y = "Sample_Sentrix_ID",all.y=T)
names(epismk)[1] <- "Sample_Sentrix_ID"

## combine fully-adjusted model covariates 
full=merge(smk, bmi, by = "Sample_Name", all.x = T)
full=merge(full,simd,by = "Sample_Name", all.x=T)
full=merge(full,ea, by = "Sample_Name", all.x= T)
full=merge(full,cells, by = "Sample_Name",all.x= T)
full=merge(full,alc,by="Sample_Name", all.x=T)
full=merge(full,epismk,by="Sample_Name", all.x=T)

## truncate to those with DNAm data - found by cell counts estimated from DNAm
full=full[-which(is.na(full$Bcell)),]

## Add in genetic PCs 
pcs=read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/GS20K_ALL_MAF5_PCA.eigenvec",header=F)
full=merge(full, pcs[,c(2,3:ncol(pcs))], by.x="Sample_Name", by.y="V2",all.x=T)

## Write out for lme/coxme models 
write.csv(full, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/covariates.csv",row.names=F)

################################
########## FACTOR ##############
################################

## Create factor covariate file for bod adjustment 
fact_cov <- data.frame(FID =samps1$Sample_Sentrix_ID,
                       IID = samps1$Sample_Sentrix_ID,
                       sex = samps1$sex,
                       Batch = samps1$Batch)
write.table(fact_cov, file="/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/bod_factors.cov", row.names=F, sep=' ', quote = T)

## Create factor covariate file for bod adjustment - no sex
fact_cov <- data.frame(FID =samps1$Sample_Sentrix_ID,
                       IID = samps1$Sample_Sentrix_ID,
                       Batch = samps1$Batch)
write.table(fact_cov, file="/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/bod_no_sex_factors.cov", row.names=F, sep=' ', quote = T)


## Create factor covariate file for EWAS 
fact_cov <- data.frame(FID =samps1$Sample_Sentrix_ID,
                       IID = samps1$Sample_Sentrix_ID,
                       sex = samps1$sex)
write.table(fact_cov, file="/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/basic_factors.cov", row.names=F, sep=' ', quote = T)

## Create factor covariate file for EWAS 
## combine with alcohol information - factor of usual 
samps2=merge(samps1,alc[,c("Sample_Name", "usual")],by="Sample_Name",all.x=T)

fact_cov <- data.frame(FID =samps2$Sample_Sentrix_ID,
                       IID = samps2$Sample_Sentrix_ID,
                       usual = samps2$usual)
write.table(fact_cov, file="/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/full_factors.cov", row.names=F, sep=' ', quote = T)



################################
########## QUANTITATIVE ########
################################

## Create quantitative covariate file for bod adjustment 
qcov <- data.frame(FID =samps1$Sample_Sentrix_ID,
                       IID = samps1$Sample_Sentrix_ID,
                       age = samps1$age)
write.table(qcov, file="/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/bod_quant.qcov", row.names=F, sep=' ', quote = T)

## Create basic quantitative covariate file for EWAS

## merge in quantitative covariates other than age 
samps2 = merge(samps1, full, by = "Sample_Name", all.x=T)

qcov <- data.frame(FID =samps2$Sample_Sentrix_ID.x,
                   IID = samps2$Sample_Sentrix_ID.x,
                   epismk = samps2$smokingScore,
                   ea = samps2$years,
                   simd = samps2$rank,
                   bmi = log(samps2$bmi),
                   units = samps2$units,
                   bcell = samps2$Bcell,
                   nk = samps2$NK,
                   gran = samps2$Gran,
                   cd4t = samps2$CD4T,
                   cd8t = samps2$CD8T)
write.table(qcov, file="/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/full_quant.qcov", row.names=F, sep=' ', quote = T)

qcov1 <- data.frame(FID =samps2$Sample_Sentrix_ID.x,
                   IID = samps2$Sample_Sentrix_ID.x,
                   epismk = samps2$smokingScore,
                   ea = samps2$years,
                   simd = samps2$rank,
                   bmi = log(samps2$bmi),
                   units = samps2$units,
                   bcell = samps2$Bcell,
                   nk = samps2$NK,
                   gran = samps2$Gran,
                   cd4t = samps2$CD4T,
                   cd8t = samps2$CD8T,
                   pc1 = samps2$V3, 
                   pc2 = samps2$V4,
                   pc3 = samps2$V5,
                   pc4 = samps2$V6,
                   pc5 = samps2$V7,
                   pc6 = samps2$V8,
                   pc7 = samps2$V9,
                   pc8 = samps2$V10,
                   pc9 = samps2$V11,
                   pc10 = samps2$V12,
                   pc11 = samps2$V13,
                   pc12 = samps2$V14,
                   pc13 = samps2$V15,
                   pc14 = samps2$V16,
                   pc15 = samps2$V17,
                   pc16 = samps2$V18,
                   pc17 = samps2$V19,
                   pc18 = samps2$V20,
                   pc19 = samps2$V21,
                   pc20 = samps2$V22)
write.table(qcov1, file="/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/pcs_quant.qcov", row.names=F, sep=' ', quote = T)


qcov2 <- data.frame(FID =samps2$Sample_Sentrix_ID.x,
                   IID = samps2$Sample_Sentrix_ID.x,
                   bcell = samps2$Bcell,
                   nk = samps2$NK,
                   gran = samps2$Gran,
                   cd4t = samps2$CD4T,
                   cd8t = samps2$CD8T)
write.table(qcov2, file="/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/wbc_quant.qcov", row.names=F, sep=' ', quote = T)







