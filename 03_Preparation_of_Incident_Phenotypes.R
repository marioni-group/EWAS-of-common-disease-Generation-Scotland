
####################################################################################

## Read in Incidence Files for each trait and combining primary and secondary cases

####################################################################################

library(gdata)
library(readxl)
library(tidyverse)
library(openxlsx)

####################################################################################

##AD 
AD <- read.csv("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/Incidence/ALZ.csv")
AD_2 <- as.data.frame(read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/AD.xlsx"))
#AD_3 <- read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/alzheimers_only_all_output_gp_with_caliber.xlsx")
AD_2 <- AD_2[,c(5,4)]
names(AD_2) <- c("id", "first")
AD_2$first <- paste0(substring(AD_2$first, 1, 4), substring(AD_2$first, 6,7))
AD <- AD[,c(1,2)]
AD = rbind(AD, AD_2)
AD = AD[order(AD$first),]
AD <- AD[which(AD$first > 199001),]
AD = AD[-which(duplicated(AD$id)),]

write.csv(AD, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/AD_combined.csv", row.names = F)

##COPD 
COPD <- read.csv("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/Incidence/COPD.csv")
COPD_2 <- as.data.frame(read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/COPD.xlsx"))
#COPD_2 <- read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/COPD_output_gp_with_caliber.xlsx")
COPD_2 = COPD_2[-grep("Acute",COPD_2$description),]
COPD_2 = COPD_2[-which(COPD_2$description %in% "Chest infection NOS"),]
COPD_2 = COPD_2[-grep("pneumonia|Pneumonia", COPD_2$description),]
COPD_2 <- COPD_2[,c(5,4)]
names(COPD_2) <- c("id", "first")
COPD_2$first <- paste0(substring(COPD_2$first, 1, 4), substring(COPD_2$first, 6,7))
COPD <- COPD[,c(1,2)]
COPD = rbind(COPD, COPD_2)
COPD = COPD[order(COPD$first),]
COPD <- COPD[which(COPD$first > 199001),]
COPD = COPD[-which(duplicated(COPD$id)),]

write.csv(COPD, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/COPD_combined.csv", row.names = F)

## Stroke 
Stroke <- read.csv("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/Incidence/Stroke.csv")
Stroke_2 <- as.data.frame(read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Stroke.xlsx"))
#Stroke_2 <- read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/stroke_output_gp_with_caliber.xlsx")
Stroke_2$start <- openxlsx::convertToDate(Stroke_2$start)
Stroke_2 <- Stroke_2[-grep("Injury|injury", Stroke_2$description),]
Stroke_2 <- Stroke_2[-grep("Mitochondrial", Stroke_2$description),]
Stroke_2 <- Stroke_2[-grep("Personal", Stroke_2$description),]
Stroke_2 <- Stroke_2[,c(5,4)]
names(Stroke_2) <- c("id", "first")
Stroke_2$first <- paste0(substring(Stroke_2$first, 1, 4), substring(Stroke_2$first, 6,7))
Stroke <- Stroke[,c(1,2)]
Stroke = rbind(Stroke, Stroke_2)
Stroke = Stroke[order(Stroke$first),]
Stroke <- Stroke[which(Stroke$first > 199001),]
Stroke = Stroke[-which(duplicated(Stroke$id)),]

write.csv(Stroke, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/Stroke_combined.csv", row.names = F)

## Lung.Cancer 
Lung.Cancer <- read.csv("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/Incidence/Lung.csv")
Lung.Cancer_2 <- as.data.frame(read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Lung_cancer.xlsx"))
#Lung.Cancer_2 <- read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/lung_cancer_output_gp_with_caliber.xlsx")
Lung.Cancer_2 <- Lung.Cancer_2[,c(5,4)]
names(Lung.Cancer_2) <- c("id", "dt")
Lung.Cancer_2$dt <- paste0(substring(Lung.Cancer_2$dt, 1, 4), substring(Lung.Cancer_2$dt, 6,7))
Lung.Cancer <- Lung.Cancer[,c(1,2,4)]
Lung.Cancer_2$smr <- 6
Lung.Cancer = rbind(Lung.Cancer, Lung.Cancer_2)
Lung.Cancer = Lung.Cancer[order(Lung.Cancer$dt),]
Lung.Cancer = Lung.Cancer[which(Lung.Cancer$dt > 199001),]
Lung.Cancer = Lung.Cancer[-which(duplicated(Lung.Cancer$id)),]

write.csv(Lung.Cancer, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/Lung_cancer_combined.csv", row.names = F)


## Breast.Cancer 
Breast.Cancer <- read.csv("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/Incidence/Breast.csv")
Breast.Cancer_2 <- as.data.frame(read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Breast_cancer.xlsx"))
#Breast.Cancer_2 <- read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/breast_cancer_output_gp_with_caliber.xlsx")
Breast.Cancer_2 <- Breast.Cancer_2[-which(Breast.Cancer_2$description %in% "Malignant neoplasm of skin of chest, excluding breast"),]
Breast.Cancer_2 <- Breast.Cancer_2[,c(5,4)]
names(Breast.Cancer_2) <- c("id", "dt")
Breast.Cancer_2$dt <- paste0(substring(Breast.Cancer_2$dt, 1, 4), substring(Breast.Cancer_2$dt, 6,7))
Breast.Cancer <- Breast.Cancer[,c(1,2,4)]
Breast.Cancer_2$smr <- 6
Breast.Cancer = rbind(Breast.Cancer, Breast.Cancer_2)
Breast.Cancer = Breast.Cancer[order(Breast.Cancer$dt),]
Breast.Cancer = Breast.Cancer[which(Breast.Cancer$dt > 199001),]
Breast.Cancer = Breast.Cancer[-which(duplicated(Breast.Cancer$id)),]
Breast.Cancer = Breast.Cancer[-which(Breast.Cancer$id %in% c(149095, 23750, 81516)),] # remove the 3 males 

write.csv(Breast.Cancer, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/Breast_cancer_combined.csv", row.names = F)

## Bowel.Cancer 
Bowel.Cancer <- read.csv("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/Incidence/Bowel.csv")
Bowel.Cancer_2 <- as.data.frame(read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Bowel_cancer.xlsx"))
#Bowel.Cancer_2 <- read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_July2020/Additional_GP_traits_extracted_4_aug_20/Bowel_cancer_GP_4_aug_20.xlsx")
Bowel.Cancer_2 <- Bowel.Cancer_2[,c(5,4)]
names(Bowel.Cancer_2) <- c("id", "dt")
Bowel.Cancer_2$dt <- paste0(substring(Bowel.Cancer_2$dt, 1, 4), substring(Bowel.Cancer_2$dt, 6,7))
Bowel.Cancer <- Bowel.Cancer[,c(1,2,4)]
Bowel.Cancer_2$smr <- 6
Bowel.Cancer = rbind(Bowel.Cancer, Bowel.Cancer_2)
Bowel.Cancer = Bowel.Cancer[order(Bowel.Cancer$dt),]
Bowel.Cancer = Bowel.Cancer[which(Bowel.Cancer$dt > 199001),]
Bowel.Cancer = Bowel.Cancer[-which(duplicated(Bowel.Cancer$id)),]

write.csv(Bowel.Cancer, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/Bowel_cancer_combined.csv", row.names = F)

##Pain 
Pain <- read.csv("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/Incidence/Back.csv")
Pain_2 <- as.data.frame(read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Back_neck_pain.xlsx"))
#Pain_2 <- read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/back_neck_pain_output_gp_with_caliber.xlsx")
Pain_2$start <- openxlsx::convertToDate(Pain_2$start)
Pain_2 <- Pain_2[,c(5,4)]
names(Pain_2) <- c("id", "first")
Pain_2$first <- paste0(substring(Pain_2$first, 1, 4), substring(Pain_2$first, 6,7))
Pain <- Pain[,c(1,2)]
Pain = rbind(Pain, Pain_2)
Pain[Pain$first=="NANA","first"] <- "189912"
Pain = Pain[order(Pain$first),]
Pain = Pain[which(Pain$first > 199001),]
Pain = Pain[-which(duplicated(Pain$id)),]

write.csv(Pain, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/Pain_combined.csv", row.names = F)


# RA
RA <- read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/extraction_SMR_4_aug_20/RA_SMR_extracted_4_aug_20.xlsx")
RA <- as.data.frame(RA)
RA <- RA[,c(8,9)]
RA <- RA[c(2,1)]
names(RA) <- c("id", "first")
RA_2 <- as.data.frame(read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/RA.xlsx"))
#RA_2 <- read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/rheumatoid_arthritis_output_gp_with_caliber.xlsx")
RA_2 <- RA_2[-grep("Polyneuropathy|Still|nodule|Juvenile|juvenile|Arthropathy", RA_2$description),]
RA_2 <- RA_2[,c(5,4)]
names(RA_2) <- c("id", "first")
RA_2$first <- paste0(substring(RA_2$first, 1, 4), substring(RA_2$first, 6,7))
RA = rbind(RA, RA_2)
RA = RA[order(RA$first),]
RA = RA[which(RA$first > 199001),]
RA = RA[-which(duplicated(RA$id)),]


write.csv(RA, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/RA_combined.csv", row.names = F)



# IBD
IBD <- read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/extraction_SMR_4_aug_20/IBD_SMR_extracted_4_aug_20.xlsx")
IBD <- as.data.frame(IBD)
IBD <- IBD[,c(8,9)]
IBD <- IBD[c(2,1)]
names(IBD) <- c("id", "first")
IBD_2 <- as.data.frame(read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/IBD.xlsx"))
#IBD_2 <- read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/inflam_bowel_disease_output_gp_with_caliber.xlsx")
IBD_2$start <- openxlsx::convertToDate(IBD_2$start)
IBD_2 <- IBD_2[-grep("Arthropathy", IBD_2$description),]
IBD_2 <- IBD_2[,c(5,4)]
names(IBD_2) <- c("id", "first")
IBD_2$first <- paste0(substring(IBD_2$first, 1, 4), substring(IBD_2$first, 6,7))
IBD = rbind(IBD, IBD_2)
IBD[IBD$first=="NANA","first"] <- "189912"
IBD = IBD[order(IBD$first),]
IBD = IBD[which(IBD$first > 199001),]
IBD = IBD[-which(duplicated(IBD$id)),]


write.csv(IBD, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/IBD_combined.csv", row.names = F)

## Heart Disease 
IHD <- read.csv("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/Incidence/Heart.csv")
IHD_2 <- as.data.frame(read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/IHD.xlsx"))
#IHD_2 <- read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_July2020/Additional_GP_traits_extracted_4_aug_20/Ischaemic_heart_disease_GP_4_aug_20.xlsx")
IHD_2 <- IHD_2[,c(5,4)]
IHD_2$start <- openxlsx::convertToDate(IHD_2$start)
names(IHD_2) <- c("id", "first")
IHD_2$first <- paste0(substring(IHD_2$first, 1, 4), substring(IHD_2$first, 6,7))
IHD <- IHD[,c(1,2)]
IHD = rbind(IHD, IHD_2)
IHD[IHD$first=="NANA","first"] <- "189912"
IHD = IHD[order(IHD$first),]
IHD = IHD[which(IHD$first > 199001),]
IHD = IHD[-which(duplicated(IHD$id)),]

write.csv(IHD, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/IHD_combined.csv", row.names = F)


## Make combined CVD phenotype for Ola/Riccardo 
# Heart Disease 
IHD <- read.csv("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/Incidence/Heart.csv")
IHD_2 <- as.data.frame(read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/IHD.xlsx"))
#IHD_2 <- read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_July2020/Additional_GP_traits_extracted_4_aug_20/Ischaemic_heart_disease_GP_4_aug_20.xlsx")
IHD_2 <- IHD_2[,c(5,4)]
IHD_2$start <- openxlsx::convertToDate(IHD_2$start)
names(IHD_2) <- c("id", "first")
IHD_2$first <- paste0(substring(IHD_2$first, 1, 4), substring(IHD_2$first, 6,7))
IHD <- IHD[,c(1,2)]
IHD = rbind(IHD, IHD_2)
IHD[IHD$first=="NANA","first"] <- "189912"
# Stroke 
Stroke <- read.csv("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/Incidence/Stroke.csv")
Stroke_2 <- as.data.frame(read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Stroke.xlsx"))
#Stroke_2 <- read_excel("U:/Datastore/IGMM/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/stroke_output_gp_with_caliber.xlsx")
Stroke_2$start <- openxlsx::convertToDate(Stroke_2$start)
Stroke_2 <- Stroke_2[-grep("Injury|injury", Stroke_2$description),]
Stroke_2 <- Stroke_2[-grep("Mitochondrial", Stroke_2$description),]
Stroke_2 <- Stroke_2[-grep("Personal", Stroke_2$description),]
Stroke_2 <- Stroke_2[,c(5,4)]
names(Stroke_2) <- c("id", "first")
Stroke_2$first <- paste0(substring(Stroke_2$first, 1, 4), substring(Stroke_2$first, 6,7))
Stroke <- Stroke[,c(1,2)]
Stroke = rbind(Stroke, Stroke_2)
# Combined 
total=rbind(IHD,Stroke)
total = total[order(total$first),]
total = total[which(total$first > 199001),]
total = total[which(!duplicated(total$id)),]
write.csv(total, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/CVD_Ola_Riccardo.csv",row.names=F)


## Osteoarthritis 
prim=read.csv("U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/new_diseases_25082022/2022-08-25_GP_readcodes_for_rob.csv")
sec=read.csv("U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/new_diseases_25082022/2022-08-25_secondary_icdcodes_for_rob.csv")
prim=prim[prim$Trait %in% "Osteoarthritis",]
sec=sec[sec$Trait %in% "Osteoarthritis",]
prim1=prim[,c("id","dt1")]
names(prim1)=c("id","first")
sec1=sec[,c("id","dt1_ym")]
names(sec1)=c("id","first")
total=rbind(prim1,sec1)
total = total[order(total$first),]
total = total[which(total$first > 199001),]
total = total[which(!duplicated(total$id)),]
write.csv(total, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/Osteoarthritis_combined.csv", row.names=F)

## Liver Fibrosis/Cirrhosis 
prim=read.csv("U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/new_diseases_25082022/2022-08-25_GP_readcodes_for_rob.csv")
sec=read.csv("U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/new_diseases_25082022/2022-08-25_secondary_icdcodes_for_rob.csv")
prim=prim[prim$Trait %in% "Liver_Cirrhosis",]
sec=sec[sec$Trait %in% "Liver_Cirrhosis",]
prim1=prim[,c("id","dt1")]
names(prim1)=c("id","first")
sec1=sec[,c("id","dt1_ym")]
names(sec1)=c("id","first")
total=rbind(prim1,sec1)
total = total[order(total$first),]
total = total[which(total$first > 199001),]
total = total[which(!duplicated(total$id)),]
write.csv(total, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/Liver_Cirrhosis_combined.csv", row.names=F)

## Chronic_Kidney_Disease 
prim=read.csv("U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/new_diseases_25082022/2022-08-25_GP_readcodes_for_rob.csv")
sec=read.csv("U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/new_diseases_25082022/2022-08-25_secondary_icdcodes_for_rob.csv")
prim=prim[prim$Trait %in% "Chronic_Kidney_Disease",]
sec=sec[sec$Trait %in% "Chronic_Kidney_Disease",]
prim1=prim[,c("id","dt1")]
names(prim1)=c("id","first")
sec1=sec[,c("id","dt1_ym")]
names(sec1)=c("id","first")
total=rbind(prim1,sec1)
total = total[order(total$first),]
total = total[which(total$first > 199001),]
total = total[which(!duplicated(total$id)),]
write.csv(total, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/CKD_combined.csv", row.names=F)

## Parkinsons
prim=read.csv("U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/new_diseases_25082022/2022-08-25_GP_readcodes_for_rob.csv")
sec=read.csv("U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/new_diseases_25082022/2022-08-25_secondary_icdcodes_for_rob.csv")
prim=prim[prim$Trait %in% "Parkinsons",]
sec=sec[sec$Trait %in% "Parkinsons",]
prim1=prim[,c("id","dt1")]
names(prim1)=c("id","first")
sec1=sec[,c("id","dt1_ym")]
names(sec1)=c("id","first")
total=rbind(prim1,sec1)
total = total[order(total$first),]
total = total[which(total$first > 199001),]
total = total[which(!duplicated(total$id)),]
write.csv(total, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/Parkinsons_combined.csv", row.names=F)


## Cervical Cancer
prim=read.csv("U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/new_diseases_25082022/2022-08-25_GP_readcodes_for_rob.csv")
sec=read.csv("U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/new_diseases_25082022/2022-08-25_secondary_icdcodes_for_rob.csv")
prim=prim[prim$Trait %in% "Cervical_Cancer",]
sec=sec[sec$Trait %in% "Cervical_Cancer",]
prim1=prim[,c("id","dt1")]
names(prim1)=c("id","first")
sec1=sec[,c("id","dt1_ym")]
names(sec1)=c("id","first")
total=rbind(prim1,sec1)
total = total[order(total$first),]
total = total[which(total$first > 199001),]
total = total[which(!duplicated(total$id)),]
write.csv(total, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/Cervical_Cancer_combined.csv", row.names=F)

## Ovarian Cancer
prim=read.csv("U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/new_diseases_25082022/2022-08-25_GP_readcodes_for_rob.csv")
sec=read.csv("U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/new_diseases_25082022/2022-08-25_secondary_icdcodes_for_rob.csv")
prim=prim[prim$Trait %in% "Ovarian_Cancer",]
sec=sec[sec$Trait %in% "Ovarian_Cancer",]
prim1=prim[,c("id","dt1")]
names(prim1)=c("id","first")
sec1=sec[,c("id","dt1_ym")]
names(sec1)=c("id","first")
total=rbind(prim1,sec1)
total = total[order(total$first),]
total = total[which(total$first > 199001),]
total = total[which(!duplicated(total$id)),]
write.csv(total, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/Ovarian_Cancer_combined.csv", row.names=F)


## Prostate Cancer
prim=read.csv("U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/new_diseases_25082022/2022-08-25_GP_readcodes_for_rob.csv")
sec=read.csv("U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/new_diseases_25082022/2022-08-25_secondary_icdcodes_for_rob.csv")
prim=prim[prim$Trait %in% "Prostate_Cancer",]
sec=sec[sec$Trait %in% "Prostate_Cancer",]
prim1=prim[,c("id","dt1")]
names(prim1)=c("id","first")
sec1=sec[,c("id","dt1_ym")]
names(sec1)=c("id","first")
total=rbind(prim1,sec1)
total = total[order(total$first),]
total = total[which(total$first > 199001),]
total = total[which(!duplicated(total$id)),]
write.csv(total, "U:/Datastore/IGMM/marioni-lab/Rob/EWAS_Disease_GS/Incidence/Prostate_Cancer_combined.csv", row.names=F)
