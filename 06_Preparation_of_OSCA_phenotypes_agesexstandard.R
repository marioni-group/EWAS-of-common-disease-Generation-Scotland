## Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/")

## Load requisite libraries 
library(survival)

## Read in covariate and consent files 
cov=read.csv("covariates.csv")
cov1=cov[,-which(colnames(cov) %in% c("units", "rank", "bmi","years", "smokingScore","pack_years"))]
cons=read.csv("2022-05-17_gp_consent.csv")


##############################################################
###### INCIDENCE PHENOTYPES - ALL INCIDENCE ##################
##############################################################

## Extract incidence phenotypes
inc=list.files("Incident_Phenotypes/", ".")
inc=inc[-which(inc%in%c("10year", "5year"))]

## Process all phenotypes except for sex-specific phenotypes
inc.sex=inc[which(inc %in% c("breast_cancer.phen","cervical_cancer.phen","ovarian_cancer.phen","prostate_cancer.phen"))]
inc=inc[-which(inc %in% c("breast_cancer.phen","cervical_cancer.phen","ovarian_cancer.phen","prostate_cancer.phen"))]

## Set up loop for processing 
for(i in inc){ 
## read in file  
tmp=read.table(paste0("Incident_Phenotypes/", i),header=T)  
## Merge in with covariate file
tmp1=merge(tmp, cov1, by = "Sample_Name")
## Subset to consenting individuals 
# tmp1=tmp1[which(tmp1$Sample_Name %in% cons$id),]
## Run coxph model to get martingale residuals  
mod=coxph(Surv(tmp1$tte, tmp1$Event) ~ tmp1$age + factor(tmp1$sex))
## Extract martingale residuals 
mod.resid=as.data.frame(residuals(mod, type = "martingale"))
names(mod.resid)[1] <- "phen"
mod.resid$FID<-tmp1$Sample_Sentrix_ID.x
mod.resid$IID<-tmp1$Sample_Sentrix_ID.x
## Tidy up dataframe
mod.resid<-mod.resid[,c("FID","IID","phen")]
## Write out dataframe 
write.table(mod.resid, paste0("Incident/Phenotypes/basic/",i),row.names=F, sep=' ', quote = T)
write.table(mod.resid, paste0("Incident/Phenotypes/wbcs/",i),row.names=F, sep=' ', quote = T)
write.table(mod.resid, paste0("Incident/Phenotypes/full/",i),row.names=F, sep=' ', quote = T)
write.table(mod.resid, paste0("Incident/Phenotypes/PCs/",i),row.names=F, sep=' ', quote = T)
print(cox.zph(mod))
}

########### SEX-SPECIFIC DISEASES ##############

## Set up loop for processing 
for(i in inc.sex){ 
  ## read in file  
  tmp=read.table(paste0("Incident_Phenotypes/", i),header=T)  
  ## Merge in with covariate file
  tmp1=merge(tmp, cov1, by = "Sample_Name")
  ## Subset to consenting individuals 
  # tmp1=tmp1[which(tmp1$Sample_Name %in% cons$id),]
  ## Run coxph model to get martingale residuals  
  mod=coxph(Surv(tmp1$tte, tmp1$Event) ~ tmp1$age)
  ## Extract martingale residuals 
  mod.resid=as.data.frame(residuals(mod, type = "martingale"))
  names(mod.resid)[1] <- "phen"
  mod.resid$FID<-tmp1$Sample_Sentrix_ID.x
  mod.resid$IID<-tmp1$Sample_Sentrix_ID.x
  ## Tidy up dataframe
  mod.resid<-mod.resid[,c("FID","IID","phen")]
  ## Write out dataframe 
  write.table(mod.resid, paste0("Incident/Phenotypes/basic/",i),row.names=F, sep=' ', quote = T)
  write.table(mod.resid, paste0("Incident/Phenotypes/wbcs/",i),row.names=F, sep=' ', quote = T)
  write.table(mod.resid, paste0("Incident/Phenotypes/full/",i),row.names=F, sep=' ', quote = T)
  write.table(mod.resid, paste0("Incident/Phenotypes/PCs/",i),row.names=F, sep=' ', quote = T)
}

##############################################################
###### PREVALENT PHENOTYPES - EXCLUDING SEX-SPECIFIC #########
##############################################################

## Extract prevalent phenotypes
inc=list.files("Prevalent_Phenotypes/", ".")
inc=inc[-which(inc %in% "asthma.phen")]

## Process all phenotypes except for breast cancer - will process separately 
inc.sex=inc[which(inc %in% c("prostate_cancer.phen","breast_cancer.phen"))]
inc=inc[-which(inc %in% c("prostate_cancer.phen","breast_cancer.phen"))]

###### BASIC MODEL #######

## Set up loop for processing 
for(i in inc){ 
  ## read in file  
  tmp=read.table(paste0("Prevalent_Phenotypes/", i),header=T)  
  ## Merge in with covariate file
  tmp1=merge(tmp, cov1, by.x = "IID", by.y = "Sample_Sentrix_ID", all.x=T)
  ## Run glm model to get residuals  
  mod=glm(factor(tmp1$phen) ~ tmp1$age + factor(tmp1$sex),family="binomial",na.action="na.exclude")
  ## Extract residuals 
  mod.resid=as.data.frame(resid(mod))
  names(mod.resid)[1] <- "phen"
  mod.resid$FID<-tmp1$FID
  mod.resid$IID<-tmp1$IID
  ## Tidy up dataframe
  mod.resid<-mod.resid[,c("FID","IID","phen")]
  ## Write out dataframe 
  write.table(mod.resid, paste0("Prevalent/Phenotypes/basic/",i),row.names=F, sep=' ', quote = T)
  write.table(mod.resid, paste0("Prevalent/Phenotypes/wbcs/",i),row.names=F, sep=' ', quote = T)
  write.table(mod.resid, paste0("Prevalent/Phenotypes/full/",i),row.names=F, sep=' ', quote = T)
  write.table(mod.resid, paste0("Prevalent/Phenotypes/PCs/",i),row.names=F, sep=' ', quote = T)
  }

########### SEX-SPECIFIC DISEASES ##############

###### BASIC MODEL #######

## Set up loop for processing 
for(i in inc.sex){ 
  ## read in file  
  tmp=read.table(paste0("Prevalent_Phenotypes/", i),header=T)
  ## Merge in with covariate file
  tmp1=merge(tmp, cov1, by.x = "IID", by.y = "Sample_Sentrix_ID", all.x=T)
  ## Run glm model to get residuals  
  mod=glm(factor(tmp1$phen) ~ tmp1$age,family="binomial",na.action="na.exclude")
  ## Extract residuals 
  mod.resid=as.data.frame(resid(mod))
  names(mod.resid)[1] <- "phen"
  mod.resid$FID<-tmp1$FID
  mod.resid$IID<-tmp1$IID
  ## Tidy up dataframe
  mod.resid<-mod.resid[,c("FID","IID","phen")]
  ## Write out dataframe 
  write.table(mod.resid, paste0("Prevalent/Phenotypes/basic/",i),row.names=F, sep=' ', quote = T)
  write.table(mod.resid, paste0("Prevalent/Phenotypes/wbcs/",i),row.names=F, sep=' ', quote = T)
  write.table(mod.resid, paste0("Prevalent/Phenotypes/full/",i),row.names=F, sep=' ', quote = T)
  write.table(mod.resid, paste0("Prevalent/Phenotypes/PCs/",i),row.names=F, sep=' ', quote = T)
  }




# ## CVD for Ola/Riccardo
# tmp=read.table("CVD.phen",header=T)
# ## Merge in with covariate file
# tmp1=merge(tmp, cov1, by = "Sample_Name")
# ## Run coxph model to get martingale residuals  
# mod=coxph(Surv(tmp1$tte, tmp1$Event) ~ tmp1$age + factor(tmp1$sex))
# ## Extract martingale residuals 
# mod.resid=as.data.frame(residuals(mod, type = "martingale"))
# names(mod.resid)[1] <- "phen"
# mod.resid$FID<-tmp1$Sample_Sentrix_ID.x
# mod.resid$IID<-tmp1$Sample_Sentrix_ID.x
# ## Tidy up dataframe
# mod.resid<-mod.resid[,c("FID","IID","phen")]
# ## Write out dataframe 
# write.table(mod.resid, "CVD_osca.phen",row.names=F, sep=' ', quote = T)
