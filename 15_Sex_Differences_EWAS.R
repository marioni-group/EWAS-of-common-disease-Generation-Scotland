## Sex-specific analyses of common disease 

## In R, make male and female text files for OSCA to only keep respective sex in analyses
R
setwd("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/")
cov=read.csv("covariates.csv")
## Subset to males and females 
males=cov[cov$sex %in% "M",]
females=cov[cov$sex %in% "F",]
## Prepare final files 
  # males
mal=males[,c("Sample_Sentrix_ID","Sample_Sentrix_ID")]
names(mal)=c("FID","IID")
write.table(mal, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/males.list", row.names=F,col.names=F,sep="\t",quote=T)
# females
femal=females[,c("Sample_Sentrix_ID","Sample_Sentrix_ID")]
names(femal)=c("FID","IID")
write.table(femal, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/females.list", row.names=F,col.names=F,sep="\t",quote=T)
# quit 
q()
n

################################################
######## PREVALENCE ############################
################################################

#### BASIC - MALES ####

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/basic/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../males.list \
--pheno $i \
--fast-linear \
--out ../../../Prevalent_Sex/males/basic/${A} \
--methylation-m
done 

#### BASIC - FEMALES ####

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/basic/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../females.list \
--pheno $i \
--fast-linear \
--out ../../../Prevalent_Sex/females/basic/${A} \
--methylation-m
done 


### WBCs - MALES ###

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/wbcs/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../males.list \
--pheno $i \
--qcovar ../../../wbc_quant.qcov \
--fast-linear \
--out ../../../Prevalent_Sex/males/wbcs/${A} \
--methylation-m
done 

### WBCs - FEMALES ###

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/wbcs/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../females.list \
--pheno $i \
--qcovar ../../../wbc_quant.qcov \
--fast-linear \
--out ../../../Prevalent_Sex/females/wbcs/${A} \
--methylation-m
done 


### FULLY-ADJUSTED - MALES ###
cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/full/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../males.list \
--pheno $i \
--qcovar ../../../full_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../../Prevalent_Sex/males/full/${A} \
--methylation-m
done 

### FULLY-ADJUSTED - FEMALES ###


cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/full/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../females.list \
--pheno $i \
--qcovar ../../../full_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../../Prevalent_Sex/females/full/${A} \
--methylation-m
done 


### FULLY-ADJUSTED + PCS - MALES ###

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/PCs/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../males.list \
--pheno $i \
--qcovar ../../../pcs_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../../Prevalent_Sex/males/PCs/${A} \
--methylation-m
done 


### FULLY-ADJUSTED + PCS - FEMALES ###

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/PCs/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../females.list \
--pheno $i \
--qcovar ../../../pcs_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../../Prevalent_Sex/females/PCs/${A} \
--methylation-m
done 



################################################
######## INCIDENCE #############################
################################################


#### BASIC - MALES ####

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/basic/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../males.list \
--pheno $i \
--fast-linear \
--out ../../../Incident_Sex/males/basic/${A} \
--methylation-m
done 

#### BASIC - FEMALES ####

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/basic/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../females.list \
--pheno $i \
--fast-linear \
--out ../../../Incident_Sex/females/basic/${A} \
--methylation-m
done 


##### WBC - MALES #####

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/wbcs/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../males.list \
--pheno $i \
--qcovar ../../../wbc_quant.qcov \
--fast-linear \
--out ../../../Incident_Sex/males/wbcs/${A} \
--methylation-m
done 

##### WBC - FEMALES #####

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/wbcs/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../females.list \
--pheno $i \
--qcovar ../../../wbc_quant.qcov \
--fast-linear \
--out ../../../Incident_Sex/females/wbcs/${A} \
--methylation-m
done 


#### FULLY-ADJUSTED - MALES ###

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/full/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../males.list \
--pheno $i \
--qcovar ../../../full_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../../Incident_Sex/males/full/${A} \
--methylation-m
done 

#### FULLY-ADJUSTED - FEMALES ###

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/full/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../females.list \
--pheno $i \
--qcovar ../../../full_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../../Incident_Sex/females/full/${A} \
--methylation-m
done 


#### FULLY-ADJUSTED + PCs - MALES #####

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/PCs/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../males.list \
--pheno $i \
--qcovar ../../../pcs_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../../Incident_Sex/males/PCs/${A} \
--methylation-m
done 


#### FULLY-ADJUSTED + PCs - FEMALES #####

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/PCs/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--keep ../../../females.list \
--pheno $i \
--qcovar ../../../pcs_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../../Incident_Sex/females/PCs/${A} \
--methylation-m
done 