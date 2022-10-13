## Adjustment for DNAm for covariates 

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/
  osca_Linux --befile osca_20k --covar bod_factors.cov --qcovar bod_quant.qcov --adj-probe --make-bod --out osca_20k

## no sex 
osca_Linux --efile osca_20k --covar bod_no_sex_factors.cov --qcovar bod_quant.qcov --adj-probe --make-bod --out osca_20k_no_sex

## Make ORM 
osca_Linux --befile osca_20k --make-orm --out osca_20k_orm

## no sex
osca_Linux --befile osca_20k_no_sex --make-orm --out osca_20k_no_sex_orm

### PREVALENT #### 
## Run basic models 

screen -S b1
cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/basic/
  
for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20k \
--pheno $i \
--fast-linear \
--out ../../Outputs/basic/${A} \
--methylation-m
done 

## sex-specific diseases 
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno breast_cancer.phen \
--fast-linear \
--out ../../Outputs/basic/breast_cancer \
--methylation-m

osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno prostate_cancer.phen \
--fast-linear \
--out ../../Outputs/basic/prostate_cancer \
--methylation-m

### WBCs ###

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/wbcs/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20k \
--pheno $i \
--qcovar ../../../wbc_quant.qcov \
--fast-linear \
--out ../../Outputs/wbcs/${A} \
--methylation-m
done 

## sex-specific diseases 
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno breast_cancer.phen \
--qcovar ../../../wbc_quant.qcov \
--fast-linear \
--out ../../Outputs/wbcs/breast_cancer \
--methylation-m

osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno prostate_cancer.phen \
--qcovar ../../../wbc_quant.qcov \
--fast-linear \
--out ../../Outputs/wbcs/prostate_cancer \
--methylation-m

### FULLY-ADJUSTED ###

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/full/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20k \
--pheno $i \
--qcovar ../../../full_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/full/${A} \
--methylation-m
done 

## sex-specific diseases 
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno breast_cancer.phen \
--qcovar ../../../full_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/full/breast_cancer \
--methylation-m

osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno prostate_cancer.phen \
--qcovar ../../../full_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/full/prostate_cancer \
--methylation-m

### FULLY-ADJUSTED + PCS ###


cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/PCs/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20k \
--pheno $i \
--qcovar ../../../pcs_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/PCs/${A} \
--methylation-m
done 

## sex-specific diseases 
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno breast_cancer.phen \
--qcovar ../../../pcs_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/PCs/breast_cancer \
--methylation-m

osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno prostate_cancer.phen \
--qcovar ../../../pcs_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/PCs/prostate_cancer \
--methylation-m

## Incident Phenotypes - cox models

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/basic/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20k \
--pheno $i \
--fast-linear \
--out ../../Outputs/basic/${A} \
--methylation-m
done 


## sex-specific diseases
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno breast_cancer.phen \
--fast-linear \
--out ../../Outputs/basic/breast_cancer \
--methylation-m

osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno prostate_cancer.phen \
--fast-linear \
--out ../../Outputs/basic/prostate_cancer \
--methylation-m

osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno ovarian_cancer.phen \
--fast-linear \
--out ../../Outputs/basic/ovarian_cancer \
--methylation-m

osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno cervical_cancer.phen \
--fast-linear \
--out ../../Outputs/basic/cervical_cancer \
--methylation-m


## Fully-adjusted models 

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/wbcs/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20k \
--pheno $i \
--qcovar ../../../wbc_quant.qcov \
--fast-linear \
--out ../../Outputs/wbcs/${A} \
--methylation-m
done 

## sex-specific diseases 
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno breast_cancer.phen \
--qcovar ../../../wbc_quant.qcov \
--fast-linear \
--out ../../Outputs/wbcs/breast_cancer \
--methylation-m

osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno prostate_cancer.phen \
--qcovar ../../../wbc_quant.qcov \
--fast-linear \
--out ../../Outputs/wbcs/prostate_cancer \
--methylation-m

osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno cervical_cancer.phen \
--qcovar ../../../wbc_quant.qcov \
--fast-linear \
--out ../../Outputs/wbcs/cervical_cancer \
--methylation-m


osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno ovarian_cancer.phen \
--qcovar ../../../wbc_quant.qcov \
--fast-linear \
--out ../../Outputs/wbcs/ovarian_cancer \
--methylation-m

#### FULLY-ADJUSTED MODELS #####

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/full/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20k \
--pheno $i \
--qcovar ../../../full_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/full/${A} \
--methylation-m
done 

## sex-specific diseases  
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno breast_cancer.phen \
--qcovar ../../../full_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/full/breast_cancer \
--methylation-m

osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno prostate_cancer.phen \
--qcovar ../../../full_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/full/prostate_cancer \
--methylation-m

osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno cervical_cancer.phen \
--qcovar ../../../full_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/full/cervical_cancer \
--methylation-m

osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno ovarian_cancer.phen \
--qcovar ../../../full_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/full/ovarian_cancer \
--methylation-m

#### FULLY-ADJUSTED MODELS + PCs #####

cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/PCs/
  
  for i in *.phen
do  
A=$( echo $i | cut -d"." -f1)
osca_Linux \
--linear \
--befile ../../../osca_20k \
--pheno $i \
--qcovar ../../../pcs_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/PCs/${A} \
--methylation-m
done 

## sex-specific diseases  
osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno breast_cancer.phen \
--qcovar ../../../pcs_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/PCs/breast_cancer \
--methylation-m

osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno prostate_cancer.phen \
--qcovar ../../../pcs_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/PCs/prostate_cancer \
--methylation-m

osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno cervical_cancer.phen \
--qcovar ../../../pcs_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/PCs/cervical_cancer \
--methylation-m

osca_Linux \
--linear \
--befile ../../../osca_20_no_sex \
--pheno ovarian_cancer.phen \
--qcovar ../../../pcs_quant.qcov \
--covar ../../../full_factors.cov \
--fast-linear \
--out ../../Outputs/PCs/ovarian_cancer \
--methylation-m



# screen -S b1
# cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/
# osca_Linux \
# --linear \
# --befile osca_20k \
# --pheno CVD_osca.phen \
# --fast-linear \
# --out CVD/CVD_basic \
# --methylation-m
# 
# cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/
#   osca_Linux \
# --linear \
# --befile osca_20k \
# --pheno CVD_osca.phen \
# --qcovar wbc_quant.qcov \
# --fast-linear \
# --out CVD/CVD_wbcs \
# --methylation-m
# 
# cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/
#   osca_Linux \
# --linear \
# --befile osca_20k \
# --pheno CVD_osca.phen \
# --qcovar full_quant.qcov \
# --covar full_factors.cov \
# --fast-linear \
# --out CVD/CVD_full \
# --methylation-m
# 
# cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/
#   osca_Linux \
# --linear \
# --befile osca_20k \
# --pheno CVD_osca.phen \
# --qcovar pcs_quant.qcov \
# --covar full_factors.cov \
# --fast-linear \
# --out CVD/CVD_pcs \
# --methylation-m

