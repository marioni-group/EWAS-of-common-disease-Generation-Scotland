## Preparation of Time-To-Event variables 

## Set working directory
setwd("U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/")

## Read in incident phenotypes 
AD=read.csv("AD_combined.csv")
COPD=read.csv("COPD_combined.csv")
Diabetes=read.csv("Diabetes_combined.csv")
IBD=read.csv("IBD_combined.csv")
IHD=read.csv("IHD_combined.csv")
Pain=read.csv("Pain_combined.csv")
RA=read.csv("RA_combined.csv")
Stroke=read.csv("Stroke_combined.csv")
Bowel.Cancer=read.csv("Bowel_cancer_combined.csv")
Breast.Cancer=read.csv("Breast_cancer_combined.csv")
Lung.Cancer=read.csv("Lung_cancer_combined.csv")
Osteoarthritis=read.csv("Osteoarthritis_combined.csv")
CKD=read.csv("CKD_combined.csv")
Liver.Cirrhosis=read.csv("Liver_Cirrhosis_combined.csv")
Parkinsons=read.csv("Parkinsons_combined.csv")
Prostate.Cancer=read.csv("Prostate_Cancer_combined.csv")
Ovarian.Cancer=read.csv("Ovarian_Cancer_combined.csv")
Cervical.Cancer=read.csv("Cervical_Cancer_combined.csv")


## Read in censoring info 
cens=read.csv("../censoring_oct2020_ckd.csv")
names(cens)[127] <-"pain"


######################################################
#### Diseases with prevalent data at baseline ########
######################################################
## Set up loop 
my.list = list(COPD,Stroke,Diabetes,Pain,RA,IHD,Osteoarthritis,Prostate.Cancer,CKD)
my.names = c("COPD_Y","stroke_Y","diabetes_Y","pain","rheum_arthritis_Y","heart_disease_Y", "osteo_arthritis_Y","prostate_cancer_Y", "CKD")

for(i in 1:length(my.list)){
tmp <- my.list[[i]]
## Remove Prevalent cases 
prev.tmp = cens[-which(cens[,my.names[[i]]]%in%1),]
tmp1 = tmp[which(tmp$id %in% prev.tmp$Sample_Name),]

## Obtain Age of Onset 
affected = prev.tmp[which(prev.tmp$Sample_Name %in% tmp1$id),] 
age_onset = tmp[,c("first", "id")]
affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
affected$Event = 1
affected$yoe = substring(affected$first, 1, 4)
affected$moe = substring(affected$first, 5,6)
affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
affected$age_event = affected$age_event1 + affected$month_event1
affected$first = NULL
affected$yoe = NULL 
affected$moe = NULL
affected$month_event1 = NULL 
affected$age_event1 = NULL

healthy = prev.tmp[-which(prev.tmp$Sample_Name %in% tmp$id),]
healthy$Event = 0
healthy$age_event = 0 
names(affected)[names(affected)=="id"] <- "Sample_Name"
cox = rbind(affected, healthy)

## Prepare tte variable 
cox$age_death = 0
cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
cox$tte = cox$age_at_event - cox$Age
cox$tte = as.numeric(cox$tte)
cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
cox$Event = as.numeric(cox$Event)
cox$tte<-as.numeric(cox$tte)
cox1<-cox[,c("Sample_Name","Sample_Sentrix_ID","Event","tte")]
cox1<-cox1[which(!is.na(cox1$tte)),]

## Extract variable name and write out file
var.name=gsub("_Y", "", my.names[[i]])
write.table(cox1,paste0("../Incidence/Phenotypes_30062022/", var.name, ".phen"),row.names=F, sep=' ')

## Print to denote completion
print(var.name)
} 


######################################################
#### Diseases without prevalent data at baseline #####
######################################################
my.list = list(IBD,Liver.Cirrhosis,Cervical.Cancer,Ovarian.Cancer)
my.names = c("ibd","liver_cirrhosis", "cervical_cancer", "ovarian_cancer")
for(i in 1:length(my.list)){ 
tmp <- my.list[[i]]
## Obtain Age of Onset 
affected = cens[which(cens$Sample_Name %in% tmp$id),] 
age_onset = tmp[,c("first", "id")]
affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
affected$Event = 1
affected$yoe = substring(affected$first, 1, 4)
affected$moe = substring(affected$first, 5,6)
affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
affected$age_event = affected$age_event1 + affected$month_event1
affected$first = NULL
affected$yoe = NULL 
affected$moe = NULL
affected$month_event1 = NULL 
affected$age_event1 = NULL

healthy = cens[-which(cens$Sample_Name %in% tmp$id),]
healthy$Event = 0
healthy$age_event = 0 
names(affected)[names(affected)=="id"] <- "Sample_Name"
cox = rbind(affected, healthy)

## Prepare tte variable 
cox$age_death = 0
cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
cox$tte = cox$age_at_event - cox$Age
cox$tte = as.numeric(cox$tte)
cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
cox$Event = as.numeric(cox$Event)
cox$tte<-as.numeric(cox$tte)
cox1<-cox[,c("Sample_Name","Sample_Sentrix_ID","Event","tte")]
cox1<-cox1[which(!is.na(cox1$tte)),]

## Extract variable name and write out file
var.name=my.names[[i]]
write.table(cox1,paste0("../Incidence/Phenotypes_30062022/", var.name, ".phen"),row.names=F, sep=' ')

## Print to denote completion
print(var.name)
} 


######################################################
#### AD processed separately due to age-of-onset #####
######################################################
tmp <- AD
## Obtain Age of Onset 
affected = cens[which(cens$Sample_Name %in% tmp$id),] 
age_onset = tmp[,c("first", "id")]
affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
affected$Event = 1
affected$yoe = substring(affected$first, 1, 4)
affected$moe = substring(affected$first, 5,6)
affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
affected$age_event = affected$age_event1 + affected$month_event1
affected$first = NULL
affected$yoe = NULL 
affected$moe = NULL
affected$month_event1 = NULL 
affected$age_event1 = NULL

healthy = cens[-which(cens$Sample_Name %in% tmp$id),]
healthy$Event = 0
healthy$age_event = 0 
names(affected)[names(affected)=="id"] <- "Sample_Name"
cox = rbind(affected, healthy)

## Prepare tte variable 
cox$age_death = 0
cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
cox$tte = cox$age_at_event - cox$Age
cox$tte = as.numeric(cox$tte)
cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
cox$Event = as.numeric(cox$Event)
cox$tte<-as.numeric(cox$tte)
cox = cox[cox$age_at_event >=65,]
cox1<-cox[,c("Sample_Name","Sample_Sentrix_ID","Event","tte")]
cox1<-cox1[which(!is.na(cox1$tte)),]

## Write out file
write.table(cox1,"../Incidence/Phenotypes_30062022/ad.phen",row.names=F, sep=' ')


######################################################
#### Parkinsons processed separately due to age-of-onset #####
######################################################
tmp <- Parkinsons
## Obtain Age of Onset 
affected = cens[which(cens$Sample_Name %in% tmp$id),] 
age_onset = tmp[,c("first", "id")]
affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
affected$Event = 1
affected$yoe = substring(affected$first, 1, 4)
affected$moe = substring(affected$first, 5,6)
affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
affected$age_event = affected$age_event1 + affected$month_event1
affected$first = NULL
affected$yoe = NULL 
affected$moe = NULL
affected$month_event1 = NULL 
affected$age_event1 = NULL

healthy = cens[-which(cens$Sample_Name %in% tmp$id),]
healthy$Event = 0
healthy$age_event = 0 
names(affected)[names(affected)=="id"] <- "Sample_Name"
cox = rbind(affected, healthy)

## Prepare tte variable 
cox$age_death = 0
cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
cox$tte = cox$age_at_event - cox$Age
cox$tte = as.numeric(cox$tte)
cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
cox$Event = as.numeric(cox$Event)
cox$tte<-as.numeric(cox$tte)
cox = cox[cox$age_at_event >=65,]
cox1<-cox[,c("Sample_Name","Sample_Sentrix_ID","Event","tte")]
cox1<-cox1[which(!is.na(cox1$tte)),]

## Write out file
write.table(cox1,"../Incidence/Phenotypes_30062022/parkinsons_age65.phen",row.names=F, sep=' ')


######################################################################
#### Some cancer phenotypes processed separately due to SMR data #####
######################################################################
my.list.cancer = list(Breast.Cancer,Lung.Cancer,Bowel.Cancer)
my.names = c("breast_cancer_Y","lung_cancer_Y","bowel_cancer_Y")

l.cancer=lapply(my.list.cancer, "[", 1:3)
smr=lapply(my.list.cancer, "[", c(1,3))
smr1=lapply(smr,subset, smr==1)

for(i in 1:length(my.list.cancer)){
  tmp <- my.list.cancer[[i]]
  smr_tmp <-  smr1[[i]]
  ## Remove Prevalent cases 
  prev.tmp = cens[-which(cens[,my.names[[i]]]%in%1),]
  prev.tmp=prev.tmp[-which(prev.tmp$Sample_Name %in% smr_tmp$id),]
  tmp1 = tmp[which(tmp$id %in% prev.tmp$Sample_Name),]
  
  ## Obtain Age of Onset 
  affected = prev.tmp[which(prev.tmp$Sample_Name %in% tmp1$id),] 
  age_onset = tmp[,c("dt", "id")]
  affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
  affected$Event = 1
  affected$yoe = substring(affected$dt, 1, 4)
  affected$moe = substring(affected$dt, 5,6)
  affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
  affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
  affected$age_event = affected$age_event1 + affected$month_event1
  affected$dt = NULL
  affected$yoe = NULL 
  affected$moe = NULL
  affected$month_event1 = NULL 
  affected$age_event1 = NULL
  
  healthy = prev.tmp[-which(prev.tmp$Sample_Name %in% tmp$id),]
  healthy$Event = 0
  healthy$age_event = 0 
  names(affected)[names(affected)=="id"] <- "Sample_Name"
  cox = rbind(affected, healthy)
  
  ## Prepare tte variable 
  cox$age_death = 0
  cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
  cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
  cox$tte = cox$age_at_event - cox$Age
  cox$tte = as.numeric(cox$tte)
  cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
  cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
  cox$Event = as.numeric(cox$Event)
  cox$tte<-as.numeric(cox$tte)
  cox1<-cox[,c("Sample_Name","Sample_Sentrix_ID","Event","tte")]
  cox1<-cox1[which(!is.na(cox1$tte)),]
  
  ## Extract variable name and write out file
  var.name=gsub("_Y", "", my.names[[i]])
  write.table(cox1,paste0("../Incidence/Phenotypes_30062022/", var.name, ".phen"),row.names=F, sep=' ')
  
  ## Print to denote completion
  print(var.name)
}


## Additional QC - sex-specific phenotypes 

## Breast Cancer
breast=read.table("../Incidence/Phenotypes_30062022/breast_cancer.phen",header=T)
fem=cens[cens$sex%in%"F",]
breast1=breast[which(breast$Sample_Name %in% fem$Sample_Name),]
write.table(breast1,"../Incidence/Phenotypes_30062022/breast_cancer.phen",row.names=F, sep=' ')
## Ovarian Cancer
ovarian=read.table("../Incidence/Phenotypes_30062022/ovarian_cancer.phen",header=T)
fem=cens[cens$sex%in%"F",]
ovarian1=ovarian[which(ovarian$Sample_Name %in% fem$Sample_Name),]
write.table(ovarian1,"../Incidence/Phenotypes_30062022/ovarian_cancer.phen",row.names=F, sep=' ')
## Cervical Cancer
cervical=read.table("../Incidence/Phenotypes_30062022/cervical_cancer.phen",header=T)
fem=cens[cens$sex%in%"F",]
cervical1=cervical[which(cervical$Sample_Name %in% fem$Sample_Name),]
write.table(cervical1,"../Incidence/Phenotypes_30062022/cervical_cancer.phen",row.names=F, sep=' ')
## Prostate Cancer
prostate=read.table("../Incidence/Phenotypes_30062022/prostate_cancer.phen",header=T)
mal=cens[cens$sex%in%"M",]
prostate1=prostate[which(prostate$Sample_Name %in% mal$Sample_Name),]
write.table(prostate1,"../Incidence/Phenotypes_30062022/prostate_cancer.phen",row.names=F, sep=' ')


# ## CVD 
# tmp=read.csv("../CVD_Ola_Riccardo.csv")
# cens$CVD=0
# cens[which(cens$heart_disease_Y %in% 1 | cens$stroke_Y %in% 1), "CVD"]<-1
# 
# prev.tmp = cens[-which(cens[,"CVD"]%in%1),]
# tmp1 = tmp[which(tmp$id %in% prev.tmp$Sample_Name),]
# 
# ## Obtain Age of Onset 
# affected = prev.tmp[which(prev.tmp$Sample_Name %in% tmp1$id),] 
# age_onset = tmp[,c("first", "id")]
# affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
# affected$Event = 1
# affected$yoe = substring(affected$first, 1, 4)
# affected$moe = substring(affected$first, 5,6)
# affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
# affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
# affected$age_event = affected$age_event1 + affected$month_event1
# affected$first = NULL
# affected$yoe = NULL 
# affected$moe = NULL
# affected$month_event1 = NULL 
# affected$age_event1 = NULL
# 
# healthy = prev.tmp[-which(prev.tmp$Sample_Name %in% tmp$id),]
# healthy$Event = 0
# healthy$age_event = 0 
# names(affected)[names(affected)=="id"] <- "Sample_Name"
# cox = rbind(affected, healthy)
# 
# ## Prepare tte variable 
# cox$age_death = 0
# cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
# cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
# cox$tte = cox$age_at_event - cox$Age
# cox$tte = as.numeric(cox$tte)
# cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
# cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
# cox$Event = as.numeric(cox$Event)
# cox$tte<-as.numeric(cox$tte)
# cox1<-cox[,c("Sample_Name","Sample_Sentrix_ID","Event","tte")]
# cox1<-cox1[which(!is.na(cox1$tte)),]
# write.table(cox1,"U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/CVD.phen",row.names=F, sep=' ')
