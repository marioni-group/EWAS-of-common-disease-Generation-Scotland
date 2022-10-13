# Load requisite libraries 
library(tidyverse)
library(readxl)

# Load covid data 
t <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/2022-10-03_rob_covid_data.csv")
names(t)[1] <- "Sample_Name"

# Read in data about whether participants received covid test 
g <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/2022-10-07_positive_tests_linked.csv")
g=g[which(!is.na(g$id)),]
g=g[order(g$test_date),]
g1=g[-which(duplicated(g$id)),]
names(g1)[1]="Sample_Name"
# Ask who is in DNAm sample information 
d1 <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
d2 <- left_join(g1, d1, by = "Sample_Name")
d2=d2[which(!is.na(d2$Sample_Sentrix_ID)),]

# Subset covid data to those with positive tests and DNAm information
t1=t[which(t$Sample_Name %in% d2$Sample_Name),]
t1=merge(t1,d2,by="Sample_Name",all.x=T)

# Get appointment year and month  
t1$yoa = substring(t1$appt, 1, 4)
t1$moa = substring(t1$appt, 5,6)
# Get even year and month 
t1$yoe = substring(t1$test_date, 1, 4)
t1$moe = substring(t1$test_date, 5,6)

t1$diff=as.numeric(t1$yoe)-as.numeric(t1$yoa)
t1$month_diff=(as.numeric(t1$moe)-as.numeric(t1$moa))/12
t1$Coviddiff=t1$diff+t1$month_diff
t1$covidAge=t1$Coviddiff+t1$age

# calculate mean difference 
t1_smr <- t1[which(t1$hosp_all  %in% "1"),] 
mean <- mean(t1_smr$Coviddiff, na.rm = T)
sd <- sd(t1_smr$Coviddiff)

# Regress covid phenotype on age and sex and save out new phenotype 
tmp=as.data.frame(resid(glm(t1$hosp_all ~ t1$covidAge + factor(t1$sex),family="binomial")))
names(tmp)[1]="phen"
tmp$FID=t1$Sample_Sentrix_ID
tmp$IID=t1$Sample_Sentrix_ID
tmp=tmp[,c("FID","IID","phen")]
write.table(tmp, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/basic/covid_hospitalisation.phen",row.names=F, sep=' ', quote = T)
write.table(tmp, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/wbcs/covid_hospitalisation.phen",row.names=F, sep=' ', quote = T)
write.table(tmp, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/full/covid_hospitalisation.phen",row.names=F, sep=' ', quote = T)
write.table(tmp, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/PCs/covid_hospitalisation.phen",row.names=F, sep=' ', quote = T)




##### EWAS models ####
# basic 
cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/wbcs/
  
  osca_Linux \
--linear \
--befile ../../../osca_20k \
--pheno covid_hospitalisation.phen \
--qcovar ../../../wbc_quant.qcov \
--fast-linear \
--out ../../Outputs/wbcs/covid_hospitalisation \
--methylation-m

# full 
cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/PCs/
  
  osca_Linux \
--linear \
--befile ../../../osca_20k \
--pheno covid_hospitalisation.phen \
--qcovar ../../../pcs_quant.qcov \
--fast-linear \
--out ../../Outputs/PCs/covid_hospitalisation \
--methylation-m
