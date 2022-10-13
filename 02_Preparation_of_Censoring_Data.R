library(survival)
library(kinship2)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)

####################################################################################

#### PREP PHENOTYPE FILE WITH AGE ALIVE AND AGE DEATH INFO 

####################################################################################

## Read in GS 20k age/sex base file and add prevalent data to it 
all <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/agemonths.csv")
names(all)[2] <- "age"
# all <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/agesex.csv")
names(all)[1] <- "Sample_Name"
dim(all) # 24088   3
prevalent <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/disease.csv")
names(prevalent)[1] <- "Sample_Name"
merged_prev <- merge(all, prevalent, by = "Sample_Name", all = T)
dim(merged_prev) # 24092   115

## PREP SURVIVAL INFO

age_dead1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_COVID/Diabetes_update_codes_yipeng/2021-09-27_death.csv")
age_dead <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_COVID/Diabetes_update_codes_yipeng/2022-03-10_age_at_death.csv")
length(which(age_dead$dod_ym > 202010)) # 227 - number of individuals who died from oct 2020 onwards

## extract year/month of death as separate variables ##
age_dead$y_dead <- as.numeric(substr(age_dead$dod_ym, 1, 4))
age_dead$m_dead <- as.numeric(substr(age_dead$dod_ym, 5, 6))

# Assign the dead individuals as a separate subset
age_dead_include <- age_dead[which(age_dead$dod_ym <= 202010),] # 1350
age_dead_exclude <- age_dead[which(age_dead$dod_ym > 202010),] # 114 - rises to 227 people who died after oct2020 will be excluded after DOB update from daniel

# Calculate a more exact estimate (by year and month) for age of death in the included 1350 individuals
age_dead_include$y_diff <- age_dead_include$y_dead - age_dead_include$yob
age_dead_include$m_diff <- (age_dead_include$m_dead - age_dead_include$mob)/12
age_dead_include$diff <- age_dead_include$y_diff + age_dead_include$m_diff

# Work out those who are alive (i.e. not in the list of 1350 dead people from above)
# age_alive <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_COVID/Diabetes_update_codes_yipeng/2022-02-16_age_alive.csv")
age_alive <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_COVID/Diabetes_update_codes_yipeng/2022-03-10_age_alive.csv")
age_alive[is.na(age_alive$mob),"mob"]<-06
age_alive <- age_alive[c(1:3,6,8)]
age_alive = age_alive[!age_alive$id %in% age_dead_include$id,] # 22738 individuals who are not in the dead people we include - this covers the remainder of GS who did not die



# Ensure that all individuals in the 22738 sample are coded as alive and all individuals in the dead file are coded as such 
age_alive$dead <- 0
age_dead_include$dead <- 1

# Find age the 'alive' people were in oct 2020
age_alive$y_diff <- 2020 - age_alive$yob
age_alive$m_diff <- (10 - age_alive$mob)/12
age_alive$diff <- age_alive$y_diff + age_alive$m_diff

# The included dead people will have their age at death taken forward (i.e. pre 2020) as the aged column - this is already provided in the file 
# The excluded dead people will be classed as alive and will have their age at 2020 taken forward as the aged column - we have just calculated this as part of the wider group of alive individuals

# Subset to just the cols needed for joining dead and alive data
age_alive <- age_alive[c(1:4,8)] 
age_dead_include <- age_dead_include[c(1,2,3,7,13)]


# Bind the rows of the alive and dead files together for the whole GS sample 
names(age_dead_include) <- c("id", "yob", "mob", "dead", "aged")
names(age_alive) <- c("id", "yob", "mob", "dead", "aged")
age = rbind(age_alive, age_dead_include)
dim(age) # 24088     5
table(age$dead)

#     0     1
# 22738  1350

names(age)[1] <- "Sample_Name"

## Add survival info to the 20k base file 
d1 <- left_join(merged_prev, age, by = "Sample_Name")
dim(d1) # [1] 24092   119

# Create a subset with DOBMOB to use to filter cases by 
d2 <- d1[c(1,2,116,117)]

names(d1)[2] <- "Age"
## Clean up dataset 
d1$Age2 <- gsub(" ", "",format(round(d1$Age, 3), nsmall = 2))
d1$age_year <- as.numeric(as.data.frame(do.call(rbind, strsplit(as.character(d1$Age2),"\\.")))$V1)
d1$age_month <- as.numeric(paste0("0.", as.data.frame(do.call(rbind, strsplit(as.character(d1$Age2),"\\.")))$V2))*12

d1$base_year <- d1$yob+as.numeric(d1$age_year)
d1$base_month <- d1$mob+as.numeric(d1$age_month)

d1$base_month <- round(d1$base_month)


d1$base_year <- ifelse(d1$base_month >= 12, d1$base_year+1, d1$base_year)
d1$base_month <- ifelse(d1$base_month >= 12, d1$base_month-12, d1$base_month)

d1$base_month <- ifelse(d1$base_month < 10, paste0("0", d1$base_month), d1$base_month)
d1$base_date <- paste0(d1$base_year, d1$base_month)

d1[d1$base_date == "NANA","base_date"] <- NA

### ADD EXTRA BASELINE VARIABLES INTO DATSAET 

# Now load in pain file 
pain=read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/chronic_painv2.csv") 
pain1=read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/chronic_painv5.csv")

## combine v2 and v5 questionnaires and remove duplicates 
pain = rbind(pain[,c("ID","pain_3_months", "back_pain", "neck_pain")],pain1[,c("ID","pain_3_months", "back_pain","neck_pain")])
pain = pain[-which(duplicated(pain$ID)),]

## clean up variables 
pain[pain$back_pain %in% 2,"back_pain"]<-1
pain[pain$neck_pain %in% 2,"neck_pain"]<-1

# Assign pain from the pain file into the pehnotype file pain variable (it doesnt exist yet so needs to be added)
# First up just give everyone a 0
d1$Pain <- 0
# Then we give everyone that has back or neck pain in the reference doc and 1 - this is based on the secondary care pain file 
d1[d1$Sample_Name %in% pain[(pain$back_pain %in% 1 | pain$neck_pain %in% 1),"ID"],"Pain"] <- 1

# Do the same for AD (it was split into males and females)
d1$AD <- 0 
d1[d1$Sample_Name %in% d1[(d1$alzheimers_M %in% 1 | d1$alzheimers_F %in% 1),"Sample_Name"],"AD"] <- 1

## QC diabetes 
## Diabetes QC - keep only type 2, remove known type 1 cases or other 
diab.qc=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Diabetes.csv")
diab.qc=diab.qc[-which(diab.qc$tname %in% "Type 2"),]
d1[which(d1$ID %in% diab.qc$id),"diabetes_Y"] <- NA

## Add in DNAm basename
targs=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
targs=targs[,c("Sample_Name","Sample_Sentrix_ID")]
d1<-merge(targs,d1,by="Sample_Name")


write.csv(d1, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/censoring_oct2020.csv",row.names =F)
