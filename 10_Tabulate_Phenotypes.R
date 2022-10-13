##### Tabulate data #####

## Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/")

####################################
####### COVARIATES #################
####################################

## Read in covariate info
cov=read.csv("covariates.csv")
## Get mean and SD for each phenotype of interest 
phenos=c("age","rank","units","smokingScore","bmi","Bcell","NK","Mono","Gran","CD4T", "CD8T", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22")  
cov=cov[,c(1,5,13,14,19,which(colnames(cov)%in% phenos))]
cov.sub=cov[complete.cases(cov),]

## Get matrix to store results 
mat=as.data.frame(matrix(nrow=length(phenos),ncol=4)) 
names(mat)<-c("Trait", "n", "mean", "sd")

## start loop 
for(i in 1:length(phenos)){ 
mat[i,1] <- as.character(phenos[i])
mat[i,2] <- length(which(!is.na(cov[,phenos[i]])))
mat[i,3] <- mean(cov[,phenos[i]],na.rm=T)
mat[i,4] <- sd(cov[,phenos[i]],na.rm = T)
}

## Add in median years and and %female 
## years 
mat[length(phenos)+1,1]<-"years"
mat[length(phenos)+1,2]<-length(which(!is.na(cov[,"years"])))
mat[length(phenos)+1,3]<-median(cov$years,na.rm=T)
mat[length(phenos)+1,4]<-IQR(cov$years, na.rm = T)
## sex 
mat[length(phenos)+2,1]<-"sex"
mat[length(phenos)+2,2]<-length(which(!is.na(cov[,"sex"])))
mat[length(phenos)+2,3]<-length(which(cov[,"sex"] %in% "F"))
mat[length(phenos)+2,4]<-length(which(cov[,"sex"] %in% "F"))/nrow(cov)

## store output 
mat.all <- mat

####################################
####### PREVALENCE #################
####################################

## extract files
phen=list.files("Prevalent_Phenotypes/", ".phen")

## set up matrix to store output 
mat.phen=as.data.frame(matrix(nrow=length(phen),ncol=5))
names(mat.phen)=c("Trait","Cases.Basic", "Controls.Basic", "Cases.Full", "Controls.Full")

## start loop 
for(i in 1:length(phen)){ 
# read in file
tmp=read.table(paste0("Prevalent_Phenotypes/", phen[i]),header=T)
# extract trait name 
mat.phen[i,1] <- gsub(".phen.*", "", phen[i])
### BASIC MODEL ####
# extract no. of cases 
mat.phen[i,2] <- length(which(tmp$phen %in% 1))
# extract no. of controls
mat.phen[i,3] <- length(which(tmp$phen %in% 0))
### FULL MODEL ####
tmp.merge=merge(cov,tmp,by.x="Sample_Sentrix_ID",by.y="FID",all.x=T)
tmp.merge=tmp.merge[complete.cases(tmp.merge),]
# extract no. of cases 
mat.phen[i,4] <- length(which(tmp.merge$phen %in% 1))
# extract no. of controls 
mat.phen[i,5] <- length(which(tmp.merge$phen %in% 0))
# print to denote completion
print(i)
}

mat.phen=mat.phen[which(!mat.phen$Trait %in% c("asthma","parkinsons_family")),]
mat.phen$Trait=c("Alzheimer's Disease**", "Colorectal Cancer", "Breast Cancer*", "Chronic Kidney Disease", "Chronic Obstructive Pulmonary Disease", "Diabetes (Type 2)", "Ischemic Heart Disease", "Lung Cancer", "Osteoarthritis", "Chronic Neck/Back Pain", "Parkinson's Disease", "Prostate Cancer*", "Rheumatoid Arthritis", "Stroke")
mat.phen=mat.phen[order(mat.phen$Trait),]
write.csv(mat.phen, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/prevalent_phenotypes_count.csv", row.names=F)
####################################
####### INCIDENCE ##################
####################################

## extract files
phen=list.files("Incident_Phenotypes/", ".phen")

## set up matrix to store output 
mat.inc=as.data.frame(matrix(nrow=length(phen),ncol=9))
names(mat.inc)=c("Trait","Cases.Basic", "Controls.Basic", "mean.tte.basic", "median.tte.basic", "Cases.Full", "Controls.Full", "mean.tte.full", "median.tte.full")

## start loop 
for(i in 1:length(phen)){ 
  # read in file
  tmp=read.table(paste0("Incident_Phenotypes/", phen[i]),header=T)
  # subset to consenting individuals only 
  # tmp=tmp[which(tmp$Sample_Name %in% cons$id),]
  # extract trait name 
  mat.inc[i,1] <- gsub(".phen.*", "", phen[i])
  ### BASIC MODEL ####
  # extract no. of cases 
  mat.inc[i,2] <- length(which(tmp$Event %in% 1))
  # extract no. of controls
  mat.inc[i,3] <- length(which(tmp$Event %in% 0))
  ## mean tte 
  mat.inc[i,4] <- paste(signif(mean(tmp[tmp$Event %in% 1,"tte"]),2), paste0("(", signif(sd(tmp[tmp$Event %in% 1,"tte"]),2), ")"))
  ## median tte 
  mat.inc[i,5] <- paste(signif(median(tmp[tmp$Event %in% 1,"tte"]),2), paste0("(", signif(IQR(tmp[tmp$Event %in% 1,"tte"]),2), ")"))
  ### FULL MODEL ####
  tmp.merge=merge(cov.sub,tmp,by="Sample_Sentrix_ID",all.x=T)
  tmp.merge=tmp.merge[complete.cases(tmp.merge),]
  # extract no. of cases 
  mat.inc[i,6] <- length(which(tmp.merge$Event %in% 1))
  # extract no. of controls 
  mat.inc[i,7] <- length(which(tmp.merge$Event %in% 0))
  ## mean tte 
  mat.inc[i,8] <- paste(signif(mean(tmp.merge[tmp.merge$Event %in% 1,"tte"]),2), paste0("(", signif(sd(tmp.merge[tmp.merge$Event %in% 1,"tte"]),2), ")"))
  ## median tte 
  mat.inc[i,9] <- paste(signif(median(tmp.merge[tmp.merge$Event %in% 1,"tte"]),2), paste0("(", signif(IQR(tmp.merge[tmp.merge$Event %in% 1,"tte"]),2), ")"))
  # print to denote completion
  print(i)
}

mat.inc=mat.inc[which(!mat.inc$Trait %in% c("cervical_cancer","parkinsons_age65")),]
mat.inc$Trait=c("Alzheimer's Disease", "Colorectal Cancer", "Breast Cancer*", "Chronic Kidney Disease", "Chronic Obstructive Pulmonary Disease", "Diabetes (Type 2)", "Ischemic Heart Disease", "Inflammatory Bowel Disease", "Liver Cirrhosis", "Lung Cancer", "Osteoarthritis", "Ovarian Cancer*", "Chronic Neck/Back Pain", "Parkinson's Disease", "Prostate Cancer*", "Rheumatoid Arthritis", "Stroke")
mat.inc=mat.inc[order(mat.inc$Trait),]
write.csv(mat.inc, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/incident_phenotypes_count.csv", row.names=F)

##############################################
####### INCIDENCE - 5 years ##################
##############################################

## extract files
phen=list.files("Incident_Phenotypes/5year/", ".phen")

## set up matrix to store output 
mat.inc.5=as.data.frame(matrix(nrow=length(phen),ncol=9))
names(mat.inc.5)=c("Trait","Cases.Basic", "Controls.Basic", "mean.tte.basic", "median.tte.basic", "Cases.Full", "Controls.Full", "mean.tte.full", "median.tte.full")

## start loop 
for(i in 1:length(phen)){ 
  # read in file
  tmp=read.table(paste0("Incident_Phenotypes/5year/", phen[i]),header=T)
  # subset to consenting individuals only 
  # tmp=tmp[which(tmp$Sample_Name %in% cons$id),]
  # extract trait name 
  mat.inc.5[i,1] <- gsub(".phen.*", "", phen[i])
  ### BASIC MODEL ####
  # extract no. of cases 
  mat.inc.5[i,2] <- length(which(tmp$Event %in% 1))
  # extract no. of controls
  mat.inc.5[i,3] <- length(which(tmp$Event %in% 0))
  ## mean tte 
  mat.inc.5[i,4] <- paste(signif(mean(tmp[tmp$Event %in% 1,"tte"]),2), paste0("(", signif(sd(tmp[tmp$Event %in% 1,"tte"]),2), ")"))
  ## median tte 
  mat.inc.5[i,5] <- paste(signif(median(tmp[tmp$Event %in% 1,"tte"]),2), paste0("(", signif(IQR(tmp[tmp$Event %in% 1,"tte"]),2), ")"))
  ### FULL MODEL ####
  tmp.merge=merge(cov.sub,tmp,by="Sample_Sentrix_ID",all.x=T)
  tmp.merge=tmp.merge[complete.cases(tmp.merge),]
  # extract no. of cases 
  mat.inc.5[i,6] <- length(which(tmp.merge$Event %in% 1))
  # extract no. of controls 
  mat.inc.5[i,7] <- length(which(tmp.merge$Event %in% 0))
  ## mean tte 
  mat.inc.5[i,8] <- paste(signif(mean(tmp.merge[tmp.merge$Event %in% 1,"tte"]),2), paste0("(", signif(sd(tmp.merge[tmp.merge$Event %in% 1,"tte"]),2), ")"))
  ## median tte 
  mat.inc.5[i,9] <- paste(signif(median(tmp.merge[tmp.merge$Event %in% 1,"tte"]),2), paste0("(", signif(IQR(tmp.merge[tmp.merge$Event %in% 1,"tte"]),2), ")"))
  # print to denote completion
  print(i)
}


mat.inc.5=mat.inc.5[which(!mat.inc.5$Trait %in% c("cervical_cancer","parkinsons_age65")),]
mat.inc.5$Trait=c("Alzheimer's Disease", "Colorectal Cancer", "Breast Cancer*", "Chronic Kidney Disease", "Chronic Obstructive Pulmonary Disease", "Diabetes (Type 2)", "Ischemic Heart Disease", "Inflammatory Bowel Disease", "Liver Cirrhosis", "Lung Cancer", "Osteoarthritis", "Ovarian Cancer*", "Chronic Neck/Back Pain", "Parkinson's Disease", "Prostate Cancer*", "Rheumatoid Arthritis", "Stroke")
mat.inc.5=mat.inc.5[order(mat.inc.5$Trait),]
write.csv(mat.inc.5, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/incident_5year_phenotypes_count.csv", row.names=F)


##############################################
####### INCIDENCE - 10 years ##################
##############################################

## extract files
phen=list.files("Incident_Phenotypes/10year/", ".phen")

## set up matrix to store output 
mat.inc.10=as.data.frame(matrix(nrow=length(phen),ncol=9))
names(mat.inc.10)=c("Trait","Cases.Basic", "Controls.Basic", "mean.tte.basic", "median.tte.basic", "Cases.Full", "Controls.Full", "mean.tte.full", "median.tte.full")

## start loop 
for(i in 1:length(phen)){ 
  # read in file
  tmp=read.table(paste0("Incident_Phenotypes/10year/", phen[i]),header=T)
  # subset to consenting individuals only 
  # tmp=tmp[which(tmp$Sample_Name %in% cons$id),]
  # extract trait name 
  mat.inc.10[i,1] <- gsub(".phen.*", "", phen[i])
  ### BASIC MODEL ####
  # extract no. of cases 
  mat.inc.10[i,2] <- length(which(tmp$Event %in% 1))
  # extract no. of controls
  mat.inc.10[i,3] <- length(which(tmp$Event %in% 0))
  ## mean tte 
  mat.inc.10[i,4] <- paste(signif(mean(tmp[tmp$Event %in% 1,"tte"]),2), paste0("(", signif(sd(tmp[tmp$Event %in% 1,"tte"]),2), ")"))
  ## median tte 
  mat.inc.10[i,5] <- paste(signif(median(tmp[tmp$Event %in% 1,"tte"]),2), paste0("(", signif(IQR(tmp[tmp$Event %in% 1,"tte"]),2), ")"))
  ### FULL MODEL ####
  tmp.merge=merge(cov.sub,tmp,by="Sample_Sentrix_ID",all.x=T)
  tmp.merge=tmp.merge[complete.cases(tmp.merge),]
  # extract no. of cases 
  mat.inc.10[i,6] <- length(which(tmp.merge$Event %in% 1))
  # extract no. of controls 
  mat.inc.10[i,7] <- length(which(tmp.merge$Event %in% 0))
  ## mean tte 
  mat.inc.10[i,8] <- paste(signif(mean(tmp.merge[tmp.merge$Event %in% 1,"tte"]),2), paste0("(", signif(sd(tmp.merge[tmp.merge$Event %in% 1,"tte"]),2), ")"))
  ## median tte 
  mat.inc.10[i,9] <- paste(signif(median(tmp.merge[tmp.merge$Event %in% 1,"tte"]),2), paste0("(", signif(IQR(tmp.merge[tmp.merge$Event %in% 1,"tte"]),2), ")"))
  # print to denote completion
  print(i)
}

mat.inc.10=mat.inc.10[which(!mat.inc.10$Trait %in% c("cervical_cancer","parkinsons_age65")),]
mat.inc.10$Trait=c("Alzheimer's Disease", "Colorectal Cancer", "Breast Cancer*", "Chronic Kidney Disease", "Chronic Obstructive Pulmonary Disease", "Diabetes (Type 2)", "Ischemic Heart Disease", "Inflammatory Bowel Disease", "Liver Cirrhosis", "Lung Cancer", "Osteoarthritis", "Ovarian Cancer*", "Chronic Neck/Back Pain", "Parkinson's Disease", "Prostate Cancer*", "Rheumatoid Arthritis", "Stroke")
mat.inc.10=mat.inc.10[order(mat.inc.10$Trait),]
write.csv(mat.inc.10, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/incident_10year_phenotypes_count.csv", row.names=F)

## Save output files 
mat.all$Trait=c("Age", "SIMD", "Alcohol Consumption", "EpiSmokEr","Body Mass Index", "B Cell", "Natural Killer", "Monocytes", "Granulocytes", "CD4T+", "CD8T+", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "Education", "Sex")
mat.all$units=c("Years", "Rank", "Units/week", "Score", "kg/m2", "Score", "Score", "Score","Score", "Score","Score","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component","Principal Component", "Years","None")
mat.all=mat.all[,c("Trait","units","n", "mean","sd")]
mat.all[,4:5]=signif(mat.all[,4:5],2)
write.csv(mat.all, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/Covariates_all.csv",row.names=F)
