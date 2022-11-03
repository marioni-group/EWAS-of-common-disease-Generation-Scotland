#########################################
## Epigenomics of Common Disease ########
#########################################

## Load requisite libraries 
library(coxme)
library(survival)
library(data.table)

# Creat function for outlier removal - 4SD
outlierID <- function(x, cut=4) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}

outlierTrim <- function(x, cut=4) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}


############################################################
##### STEP 1. READ IN METHYLATION FILE  ####################
############################################################

meth=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/DNAm_18413_resid_disease.csv"))
cg=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/DNAm_cpg_list.csv",header=F))
row.names(meth)<-cg$V1

meth=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/mvals.rds")

# Match order of IDs in methylation and samps file 
samps=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
ids=samps$Sample_Sentrix_ID
meth=meth[,which(colnames(meth)%in%samps$Sample_Sentrix_ID)]
meth=meth[,match(ids,colnames(meth))]

############################################################
##### STEP 2. READ IN COVARIATE FILE #######################
############################################################

# Read in covariates  file
cov<- read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/covariates.csv")

# Match order of IDs with other files 
cov=cov[match(ids,cov$Sample_Sentrix_ID),]

# Check order of IDs match with other files 
table(cov$Sample_Sentrix_ID==samps$Sample_Sentrix_ID)
row.names(meth)=ids


############################################################
##### STEP 3. READ IN COXPH VIOLATIONS #####################
############################################################

#### We will focus on fully-adjusted models #####
res=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Sensitivity/coxph/full_model_violations.csv")
meth=meth[,which(colnames(meth)%in%res$CpG)]

##############################################################
##### STEP 4. LOOP FOR COX.ZPH MODELS ########################
##############################################################

# Set up loop to loop through results files 
# Set up matrix to store results
mat.res=as.data.frame(matrix(nrow=nrow(res)*14,ncol=12))
names(mat.res) <- c("CpG", "Trait", "Years_to_Follow", "HR", "LCI", "HCI", "P", "n.cases", "n.controls", "n.total", "Cox.Zph_Local", "Cox.Zph_Global")
# Extract model name 
for(j in 1:nrow(res)){ 
  # Loop through each association in the file  
  # Extract CpG 
  cpg=res[j,"CpG"]  
  # Extract Trait 
  trait=res[j,"Trait"]

  ## Loop through each year of follow-up (1-14 years)
  for(hr in 1:14){ 
  # Read in phenotype 
  phen=read.table(paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/CoxPH_Violations/", trait, "/year_", hr, ".phen"),header=T)
  ids=samps$Sample_Sentrix_ID
  # Tidy up dataframes so that they all match and have complete data 
  phen=phen[match(ids,phen$Sample_Sentrix_ID),]
  phen=phen[complete.cases(phen),]
  cov2=merge(phen,cov,by="Sample_Sentrix_ID")
  meth1=meth[which(row.names(meth) %in% cov2$Sample_Sentrix_ID),]
  ids3=cov2$Sample_Sentrix_ID
  meth1=meth1[match(ids3,row.names(meth1)),]
  # Run model
  mod = coxph(Surv(cov2$tte, cov2$Event) ~ scale(meth1[,as.character(cpg)]) + cov2$age + factor(cov2$sex)+ cov2$Bcell + cov2$NK + cov2$Gran + cov2$CD4T + cov2$CD8T + log(cov2$bmi)  + cov2$units + cov2$years + cov2$smokingScore  + factor(cov2$usual) + cov$rank + cov2$V3 + cov2$V4 + cov2$V5 + cov2$V6 + cov2$V7 + cov2$V8 + cov2$V9 + cov2$V10 + cov2$V11 + cov2$V12 + cov2$V13 + cov2$V14 + cov2$V15 + cov2$V16 + cov2$V17 + cov2$V18 + cov2$V19 + cov2$V20 + cov2$V21 + cov2$V22)
  # Store CpG 
  mat.res[(j+(13*(j-1)))+(hr-1),1] <- as.character(cpg)
  # Store Trait
  mat.res[(j+(13*(j-1)))+(hr-1),2] <- as.character(trait)
  # Store time to follow-up
  mat.res[(j+(13*(j-1)))+(hr-1),3] <- as.character(hr)
  # Store HR 
  mat.res[(j+(13*(j-1)))+(hr-1),4] <- signif(summary(mod)$coefficients[1,2],2)
  # Store HR 2.5 CI and 9.5 CI
  mat.res[(j+(13*(j-1)))+(hr-1),5:6] <- exp(confint(mod))[1,c(1,2)]
  # Store P 
  mat.res[(j+(13*(j-1)))+(hr-1),7] <- signif(summary(mod)$coefficients[1,5],2)
  # Store no. of cases 
  mat.res[(j+(13*(j-1)))+(hr-1),8] <- mod$nevent
  # Store no. of controls 
  mat.res[(j+(13*(j-1)))+(hr-1),9] <- mod$n-mod$nevent
  # Store total n 
  mat.res[(j+(13*(j-1)))+(hr-1),10] <- mod$n
  # Store Local P 
  mat.res[(j+(13*(j-1)))+(hr-1),11] <- as.numeric(cox.zph(mod)$table[,"p"])[1]
  # Store Global P 
  mat.res[(j+(13*(j-1)))+(hr-1),12] <- as.numeric(cox.zph(mod)$table[,"p"][nrow(cox.zph(mod)$table)])
  # Print to denote completion
  }
} 
write.csv(mat.res,"/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Sensitivity/iterate/incident_full.csv",row.names=F) 


## Outliers - cpg 2 
meth1[which(row.names(meth1) %in% phen[phen$Event==1,][which.min(cox.zph(mod)$y),"Sample_Sentrix_ID"]),cpg]=NA
