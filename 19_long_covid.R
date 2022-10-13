# Load requisite libraries 
library(tidyverse)
library(readxl)

# Load covid data
t <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/GS_CL3_Samples_LongCovid.csv")

# Filter by removing the NAs for symptom indications
d <- t %>% filter(t$S3_SymptomLength1st != "NA") # 338 remaining 

# Construct the binary variable of < 4 and > 4 weeks 
d$binary <- ifelse(d$S3_SymptomLengthAll == 1 | d$S3_SymptomLengthAll == 2, 0, 1)
names(d)[1]="Sample_Name"

# Read in DNAm sample information 
d1 <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
d2 <- left_join(d, d1, by = "Sample_Name")

#################################

### Add in the covidage calculation 
g <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/C19_test_dates_25Oct2021.csv")
names(g)[1] <- "Sample_Name"
e <- d2
e2 <- merge(e, g, by = 'Sample_Name', all.x = TRUE)

hasAntiBD <- !is.na(e2$S3_AntiBD_Year)
hasSwab <- !is.na(e2$S3_SwabDate_Year)
hasLinkedTest <- !is.na(e2$Linked_TestDate)

# Remove individuals who reported having covid in CL1 but not in CL2
Mismatch <- (e2$Had_COVID > 0) & (!e2$S2_Had_COVID > 0)

# Replace NA with 0 in covidLife1And2Mismatch
Mismatch[is.na(Mismatch)] <- 0
e2 <- e2[!Mismatch,] 
table(e2$binary)

# Read in appt table to extract all baseline appointment dates (not just those in CovidLife)
apptTable <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/2021-07-30_appt.csv")

apptToTargetTableIndex <- match(e2$Sample_Name, apptTable$id)
apptDateString <- apptTable[apptToTargetTableIndex, 'appt']

# Extract Date from first half of the timestamp
apptDate <- lapply(apptDateString, function(x) {as.Date(strsplit(x, ' ')[[1]][[1]], '%Y-%m-%d')})

# Use date from the following sources if available: linked test, antibd, swab. Else use an approximate date of 01/01/2021
covidDate <- lapply(1:nrow(e2), function(rowName) {
  row <- e2[rowName, ]
  if (!is.na(row$Linked_TestDate)) {
    as.Date(row$Linked_TestDate, '%Y-%m-%d')
  } else if (!is.na(row$S3_AntiBD_Year)) {
    as.Date(paste(row$S3_AntiBD_Year, row$S3_AntiBD_Month, row$S3_AntiBD_Day, sep = '/'), '%Y/%m/%d')
  } else if (!is.na(row$S3_SwabDate_Year)) {
    as.Date(paste(row$S3_SwabDate_Year, row$S3_SwabDate_Month, row$S3_SwabDate_Day, sep = '/'), '%Y/%m/%d')
  } else {
    as.Date('2021-01-01', '%Y-%m-%d')
  }
})

# Difference between appointment (baseline) date and covid date
covidApptDiff <- sapply(1:length(apptDate), function(i) {as.numeric(covidDate[[i]] - apptDate[[i]]) / 365})

# Add covid appt difference onto baseline age
e2$covidAge <- e2$age + covidApptDiff

# Assign difference as extra column
e2$covidDiff <- covidApptDiff

# Remove individuals without methylation 
e2=e2[which(!is.na(e2$Sample_Sentrix_ID)),]

# filter to just smr cases 
e2_smr <- e2[which(e2$binary %in% "1"),] # 56 cases

# calculate mean difference 
mean <- mean(e2_smr$covidDiff, na.rm = T) # 11.2088
sd <- sd(e2_smr$covidDiff) # 1.2365

# Regress covid phenotype on age and sex and save out new phenotype 
tmp=as.data.frame(resid(glm(e2$binary ~ e2$covidAge + factor(e2$sex),family="binomial")))
names(tmp)[1]="phen"
tmp$FID=e2$Sample_Sentrix_ID
tmp$IID=e2$Sample_Sentrix_ID
tmp=tmp[,c("FID","IID","phen")]
write.table(tmp, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/basic/long_covid.phen",row.names=F, sep=' ', quote = T)
write.table(tmp, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/wbcs/long_covid.phen",row.names=F, sep=' ', quote = T)
write.table(tmp, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/full/long_covid.phen",row.names=F, sep=' ', quote = T)
write.table(tmp, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/PCs/long_covid.phen",row.names=F, sep=' ', quote = T)


##### EWAS models ####
# basic 
cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/wbcs/
  
osca_Linux \
--linear \
--befile ../../../osca_20k \
--pheno long_covid.phen \
--qcovar ../../../wbc_quant.qcov \
--fast-linear \
--out ../../Outputs/wbcs/long_covid \
--methylation-m
 

# full 
cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident/Phenotypes/PCs/
  
  osca_Linux \
--linear \
--befile ../../../osca_20k \
--pheno long_covid.phen \
--qcovar ../../../pcs_quant.qcov \
--fast-linear \
--out ../../Outputs/PCs/long_covid \
--methylation-m
 