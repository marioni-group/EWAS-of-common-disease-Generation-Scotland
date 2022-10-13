#########################################
## Epigenomics of Common Disease ########
#########################################

## Load requisite libraries 
library(kinship2)
library(coxme)
library(data.table)
library(limma)

## Define requisite functions
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

###################################################################
##### STEP 1. PREPARE KINSHIP MATRIX FOR RELATEDNESS ANALYSES #####
###################################################################

# Read in pedigree data 
ped = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/pedigree.csv")

# Tidy up dataframe 
ped$father <- as.numeric(ped$father)
ped$mother <- as.numeric(ped$mother)
ped$father[ped$father==0] <- NA
ped$mother[ped$mother==0] <- NA
ped$sex <- as.numeric(as.factor(ped$sex))
ped$sex[ped$sex==2] <- 0
ped$sex <- ped$sex+1

# Prepare matrix 
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))
kin_model <- kinship(kin) 


############################################################
##### STEP 2. CREATE FUNCTION TO EXTRACT MODEL RESULTS #####
############################################################

# Create function 
extract_coxme_table <- function (mod){
beta <- fixef(mod)
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}


############################################################
##### STEP 3. RESIDUALISE METHYLATION FILE AS IN OSCA ######
############################################################

# Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/Chromosomes/")
# Read in DNAm sample information
samps=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")

# # Read in methylation file 
# # by chromosome 
# for(i in 1:22){ 
#   meth=readRDS(paste0("GS20k_chr", i, "_mvals.rds"))
#   # Subset to individuals passing QC 
#   meth=meth[,which(colnames(meth) %in% samps$Sample_Sentrix_ID)]
#   # Subset to probes passing QC 
#   probes=read.table("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt", header=F)
#   meth=meth[which(row.names(meth) %in% probes$V1),]
#   # Match order of IDs in phenotype and methylation file 
#   ids=samps$Sample_Sentrix_ID
#   meth=meth[,match(ids,colnames(meth))]
#   # Check order of IDs match between methylation files 
#   table(colnames(meth)==samps$Sample_Sentrix_ID)
# 
#   # Prepare covariate matrix for regressions 
#   
#   ## Regression step - residualise for age, sex and batch 
#   design.resid <- model.matrix(~sex + age + as.factor(Batch), data=samps)
#   fit.resid <- limma::lmFit(meth, design.resid)
#   gc()
#   meth <- limma::residuals.MArrayLM(fit.resid, meth)
#   meth <- meth[!is.infinite(rowSums(meth)), ]
#   rm(fit.resid)
#   gc()
#   
#   ## Write out CpGs 
#   cpgs=as.data.frame(row.names(meth))
#   names(cpgs)[1]="CpG"
#   fwrite(cpgs, paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/CpGs/GS20k_chr", i, "_cpgs.txt"),row.names=F)
#   
#   # Save out residualised file 
#   fwrite(meth, paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/DNAm/GS20k_chr", i, "_resid_mvals.txt"),row.names=F)
#   
#   ## Remove methylation object and clean up environment 
#   rm(meth)
#   gc()
# } 

###########################################
#### STEP 4. COMBINE METHYLATION FILES ####
###########################################
# 
# Change working directory
library(data.table)
setwd("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/DNAm/")
# Extract files
files=list.files(".", ".txt")
files=files[order(files)]
# Read in files
meth <- rbindlist(lapply(files,fread))
## Write out final file
fwrite(x =as.matrix(meth), "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/DNAm_18413_resid_disease.csv", sep = ",", row.names = F, col.names = F, quote = F)

# Extract CpGs - will need this for processing final results files as the row.names get lost in methylation file in BayesR
# Extract CpGs
cpgs=list.files("../CpGs/", ".txt")
cpgs=cpgs[order(cpgs)]
# Ensure that order is same as methylation file just created
# methylation
ids=gsub("GS20k_chr", "", files)
ids=gsub("_.*", "", ids)
# cpgs
ids1=gsub("GS20k_chr", "", cpgs)
ids1=gsub("_.*", "", ids1)
# is order same?
table(ids==ids1)
# Read in cpg lists
# Change working directory
setwd("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/CpGs/")
cg <- rbindlist(lapply(cpgs,fread))
fwrite(x =as.matrix(cg), "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/DNAm_cpg_list.csv", sep = ",", row.names = F, col.names = F, quote = F)
row.names(meth)<-cg$CpG

## Read in files once created
library(data.table)
meth=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/DNAm_18413_resid_disease.csv"))
cg=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/DNAm_cpg_list.csv"))
row.names(meth)<-cg$CpG

###########################################
##### STEP 5. READ IN COVARIATE FILE ######
###########################################

# Read in covariates
cov<- read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/covariates.csv")

# Match order of IDs with other files 
ids=samps$Sample_Sentrix_ID
cov=cov[match(ids,cov$Sample_Sentrix_ID),]

# Check order of IDs match with other files 
table(cov$Sample_Sentrix_ID==colnames(meth))
table(cov$Sample_Sentrix_ID==samps$Sample_Sentrix_ID)
meth=t(meth)


###########################################
##### STEP 6. READ IN RESULTS FILES  ######
###########################################

#### We will focus on fully-adjusted models #####

# Extract fully-adjusted models 
res=list.files("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Results/Results_Final/",".csv")
res=res[grep("full",res)]

# Separate into prevalent and incident result files for reading in correct phenotype in loop
prev=res[grep("prevalent", res)]
inc=res[-grep("prevalent", res)]
inc.5=inc[grep("5years", inc)]
inc.10=inc[grep("10years", inc)]
inc1=inc[-which(inc %in% c(inc.5,inc.10))]

# Set up loop to loop through results files 
for(i in res){ 
# Read in fully-adjusted model 
tmp.res=read.csv(paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Results/Results_Final/",i))
# Set up matrix to store lmekin result 
mat.res=as.data.frame(matrix(nrow=nrow(tmp.res),ncol=8))
names(mat.res) <- c("CpG", "Trait", "Beta", "SE", "Z", "P", "OSCA_Beta", "OSCA_P")
# Extract model name 
model=as.character(gsub(".csv", "", i))

# Loop through each association in the file  
for(j in 1:nrow(tmp.res)){ 
# Extract CpG 
cpg=tmp.res[j,"Probe"]  
# Extract Trait 
trait=tmp.res[j,"trait"]

# If the result file is prevalent, run lmekin model 
if(i %in% prev){ 
# Read in phenotype
phen=read.table(paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent/Phenotypes/full/", trait, ".phen"),header=T)
ids=row.names(meth)
phen=phen[match(ids,phen$FID),]
# Run lmekin model 
mod = lmekin(phen[,"phen"] ~ meth[,as.character(cpg)] + cov$Bcell + cov$NK + cov$Gran + cov$CD4T + cov$CD8T + log(cov$bmi) + cov$rank + cov$years + cov$smokingScore + cov$pack_years + factor(cov$usual) + (1|samps$Sample_Name), varlist = kin_model*2)
mod1 <- as.data.frame(extract_coxme_table(mod))
# Store CpG 
mat.res[j,1] <- as.character(cpg)
# Store Trait
mat.res[j,2] <- as.character(trait)
# Store Beta, SE and Z  
mat.res[j,3:5] <- round(mod1[2,1:3],2)
# Store P 
mat.res[j,6] <- mod1[2,4]
## If the result file is incident (all years), run coxme model 
} else if(i %in% inc1) {
# Read in phenotype 
phen=read.table(paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/full/", trait, ".phen"),header=T)
ids=row.names(meth)
phen=phen[match(ids,phen$FID),]
# Run model
mod = lmekin(phen[,"phen"] ~ meth[,as.character(cpg)] + cov$Bcell + cov$NK + cov$Gran + cov$CD4T + cov$CD8T + log(cov$bmi) + cov$rank + cov$years + cov$smokingScore + cov$pack_years + factor(cov$usual) + (1|samps$Sample_Name), varlist = kin_model*2)
# Extract coefficients 
mod1 <- as.data.frame(extract_coxme_table(mod))
# Store CpG 
mat.res[j,1] <- as.character(cpg)
# Store Trait
mat.res[j,2] <- as.character(trait)
# Store Beta, SE and Z  
mat.res[j,3:5] <- round(mod1[2,1:3],2)
# Store P 
mat.res[j,6] <- mod1[2,4]
} ## If it is incident 5 years file, then read in appropriate file, otherwise do incidence 10 years at next step  
else if(i %in% inc.5){  
 # Read in phenotype 
phen=read.table(paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident_5years/Phenotypes/full/", trait, ".phen"),header=T)
ids=row.names(meth)
phen=phen[match(ids,phen$FID),]
# Run model 
mod = lmekin(phen[,"phen"] ~ meth[,as.character(cpg)] + cov$Bcell + cov$NK + cov$Gran + cov$CD4T + cov$CD8T + log(cov$bmi) + cov$rank + cov$years + cov$smokingScore + cov$pack_years + factor(cov$usual) + (1|samps$Sample_Name), varlist = kin_model*2)
# Extract coefficients 
mod1 <- as.data.frame(extract_coxme_table(mod))
# Store CpG 
mat.res[j,1] <- as.character(cpg)
# Store Trait
mat.res[j,2] <- as.character(trait)
# Store Beta, SE and Z  
mat.res[j,3:5] <- round(mod1[2,1:3],2)
# Store P 
mat.res[j,6] <- mod1[2,4]
} else { 
# Read in phenotype 
phen=read.table(paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident_10years/Phenotypes/full/", trait, ".phen"),header=T)
ids=row.names(meth)
phen=phen[match(ids,phen$FID),]
# Run model 
mod = lmekin(phen[,"phen"] ~ meth[,as.character(cpg)] + cov$Bcell + cov$NK + cov$Gran + cov$CD4T + cov$CD8T + log(cov$bmi) + cov$rank + cov$years + cov$smokingScore + cov$pack_years + factor(cov$usual) + (1|samps$Sample_Name), varlist = kin_model*2, na.action="na.exclude")
# Extract coefficients 
mod1 <- as.data.frame(extract_coxme_table(mod))
# Store CpG 
mat.res[j,1] <- as.character(cpg)
# Store Trait
mat.res[j,2] <- as.character(trait)
# Store Beta, SE and Z  
mat.res[j,3:5] <- round(mod1[2,1:3],2)
# Store P 
mat.res[j,6] <- mod1[2,4]
} 

# Store OSCA Beta and P results 
mat.res$OSCA_Beta <- tmp.res$b
mat.res$OSCA_P <- tmp.res$p
}
# Write out result 
write.csv(mat.res, paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Sensitivity/",model,".csv"), row.names = F)
}

