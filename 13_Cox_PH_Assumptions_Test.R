#########################################
## Epigenomics of Common Disease ########
#########################################

## Load requisite libraries 
library(coxme)
library(survival)
library(data.table)

############################################################
##### STEP 1. READ IN METHYLATION FILE  ####################
############################################################


## Read in files once created
# Both Sexes 
library(data.table)
meth=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/DNAm_18413_resid_disease.csv"))
cg=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/DNAm_cpg_list.csv",header=F))
row.names(meth)<-cg$V1

# Sex-specific analyses  
library(data.table)
meth=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/DNAm_18413_resid_nosex_disease.csv"))
cg=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/DNAm_cpg_nosex_list.csv",header=F))
row.names(meth)<-cg$V1
meth=t(meth)

# Match order of IDs in methylation and samps file 
samps=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
ids=samps$Sample_Sentrix_ID

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
##### STEP 3. READ IN RESULTS FILES ########################
############################################################

#### We will focus on fully-adjusted models #####

res=list.files("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/",".csv")
res=res[grep("full",res)]
res=res[which(res %in% c("prevalent_full.csv", "incident_full.csv"))]

##############################################################
##### STEP 4. LOOP FOR COX.ZPH MODELS ########################
##############################################################

## Both sexes 

# Set up loop to loop through results files 
# Read in fully-adjusted model 
  i=res[1]
  tmp.res=read.csv(paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/",i))
  tmp.res=tmp.res[!tmp.res$trait %in% c("cervical_cancer","covid_hospitalisation", "breast_cancer", "prostate_cancer", "ovarian_cancer"),]
# Set up matrix to store cox.zph result 
  mat.res=as.data.frame(matrix(nrow=nrow(tmp.res),ncol=6))
  names(mat.res) <- c("CpG", "Trait", "Log.Odds","P", "Cox.Zph_Local", "Cox.Zph_Global")
  # Extract model name 
  model=as.character(gsub(".csv", "", i))
for(j in 1:nrow(tmp.res)){ 
# Loop through each association in the file  
# Extract CpG 
cpg=tmp.res[j,"Probe"]  
# Extract Trait 
trait=tmp.res[j,"trait"]
  # Read in phenotype 
  phen=read.table(paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident_Phenotypes/", trait, ".phen"),header=T)
  ids=samps$Sample_Sentrix_ID
  phen=phen[match(ids,phen$Sample_Sentrix_ID),]
  phen=phen[complete.cases(phen),]
  cov2=merge(phen,cov,by="Sample_Sentrix_ID")
  meth1=meth[which(row.names(meth) %in% cov2$Sample_Sentrix_ID),]
  ids3=cov2$Sample_Sentrix_ID
  meth1=meth1[match(ids3,row.names(meth1)),]
  # Run model
  mod = coxph(Surv(cov2$tte, cov2$Event) ~ scale(meth1[,as.character(cpg)])+ cov2$age + factor(cov2$sex) + cov2$Bcell + cov2$NK + cov2$Gran + cov2$CD4T + cov2$CD8T + log(cov2$bmi)  + cov2$units + cov2$years + cov2$smokingScore  + factor(cov2$usual) + cov2$V3 + cov2$V4 + cov2$V5 + cov2$V6 + cov2$V7 + cov2$V8 + cov2$V9 + cov2$V10 + cov2$V11 + cov2$V12 + cov2$V13 + cov2$V14 + cov2$V15 + cov2$V16 + cov2$V17 + cov2$V18 + cov2$V19 + cov2$V20 + cov2$V21 + cov2$V22)
  # Store CpG 
  mat.res[j,1] <- as.character(cpg)
  # Store Trait
  mat.res[j,2] <- as.character(trait)
  # Store Beta, SE and Z  
  mat.res[j,3:4] <- signif(summary(mod)$coefficients[1,c(2,5)],2)
  # Store P 
  mat.res[j,5] <- as.numeric(cox.zph(mod)$table[,"p"][1])
  # Store Global P 
  mat.res[j,6] <- as.numeric(cox.zph(mod)$table[,"p"][nrow(cox.zph(mod)$table)])
  # Print to denote completion
  print(j)
  }
 write.csv(mat.res,"/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Sensitivity/coxph/incident_full_withsex.csv",row.names=F) 
  
  ## Sex-specific analyses
  # Set up loop to loop through results files 
  # Read in fully-adjusted model 
 tmp.res=read.csv(paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/",i))
 tmp.res=tmp.res[tmp.res$trait %in% c("breast_cancer", "prostate_cancer", "ovarian_cancer", "covid_hospitalisation"),]
 # Set up matrix to store lmekin result 
 mat.res=as.data.frame(matrix(nrow=nrow(tmp.res),ncol=6))
 names(mat.res) <- c("CpG", "Trait", "Log.Odds","P", "Cox.Zph_Local", "Cox.Zph_Global")
 # Extract model name 
 model=as.character(gsub(".csv", "", i))
 for(j in 1:nrow(tmp.res)){ 
   # Loop through each association in the file  
   # Extract CpG 
   cpg=tmp.res[j,"Probe"]  
   # Extract Trait 
   trait=tmp.res[j,"trait"]
   # Read in phenotype 
   phen=read.table(paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident_Phenotypes/", trait, ".phen"),header=T)
   ids=samps$Sample_Sentrix_ID
   phen=phen[match(ids,phen$Sample_Sentrix_ID),]
   phen=phen[complete.cases(phen),]
   cov2=merge(phen,cov,by="Sample_Sentrix_ID")
   meth1=meth[which(row.names(meth) %in% cov2$Sample_Sentrix_ID),]
   ids3=cov2$Sample_Sentrix_ID
   meth1=meth1[match(ids3,row.names(meth1)),]
   # Run model
   mod = coxph(Surv(cov2$tte, cov2$Event) ~ scale(meth1[,as.character(cpg)])+ cov2$age + cov2$Bcell + cov2$NK + cov2$Gran + cov2$CD4T + cov2$CD8T + log(cov2$bmi)  + cov2$units + cov2$years + cov2$smokingScore  + factor(cov2$usual) + cov2$V3 + cov2$V4 + cov2$V5 + cov2$V6 + cov2$V7 + cov2$V8 + cov2$V9 + cov2$V10 + cov2$V11 + cov2$V12 + cov2$V13 + cov2$V14 + cov2$V15 + cov2$V16 + cov2$V17 + cov2$V18 + cov2$V19 + cov2$V20 + cov2$V21 + cov2$V22)
   # Store CpG 
   mat.res[j,1] <- as.character(cpg)
   # Store Trait
   mat.res[j,2] <- as.character(trait)
   # Store Beta, SE and Z  
   mat.res[j,3:4] <- signif(summary(mod)$coefficients[1,c(2,5)],2)
   # Store P 
   mat.res[j,5] <- as.numeric(cox.zph(mod)$table[,"p"][1])
   # Store Global P 
   mat.res[j,6] <- as.numeric(cox.zph(mod)$table[,"p"][nrow(cox.zph(mod)$table)])
   # Print to denote completion
   print(j)
 }
 write.csv(mat.res,"/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Sensitivity/coxph/incident_full_nosex.csv",row.names=F) 
 
