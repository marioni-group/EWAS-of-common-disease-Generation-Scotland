setwd("U://Datastore//IGMM//marioni-lab//Rob//EWAS_Disease_GS//") 

## Censoring Date - all years ###

d1 <- read.csv("censoring_oct2020.csv", check.names = F)
creat <- read.csv("../Epigenetic Clocks WHO/Corr_Creatinine.csv")
creat = creat[,c(1,5)]
d1 = merge(d1, creat, by.x = "Sample_Name", by.y = "id", all.x = T)
names(d1)[129] <- "Creat"

females = d1[which(d1$sex %in% "F"),]
males = d1[which(d1$sex %in% "M"),] 

females_low <- females[which(females$Creat <= 62),] 
females_high <- females[which(females$Creat > 62),] 
females_na <- females[which(females$Creat %in% NA),] 

males_low <- males[which(males$Creat <= 80),] 
males_high <- males[which(males$Creat > 80),] 
males_na <- males[which(males$Creat %in% NA),] 

females_low$eGFR <- 141*((females_low$Creat/61.9)^-0.329)*(0.993^females_low$Age)*1.018
females_high$eGFR <- 141*((females_high$Creat/61.9)^-1.209)*(0.993^females_high$Age)*1.018
males_low$eGFR <- 141*((males_low$Creat/79.6)^-0.411)*(0.993^males_low$Age)
males_high$eGFR <-141*((males_high$Creat/79.6)^-1.209)*(0.993^males_high$Age)
females_na$eGFR <- NA
males_na$eGFR <- NA 
females_1 = rbind(females_low,females_high)
females_1 = rbind(females_1, females_na)
males_1 = rbind(males_low,males_high)
males_1 = rbind(males_1, males_na)
all = rbind(males_1, females_1)
d2<-all
d2$CKD = ifelse(d2$eGFR < 60, 1, 0)

write.csv(d2, "censoring_oct2020_ckd.csv",row.names=F)

## Censoring Date - 10 year cut-off ###

d1 <- read.csv("censoring_10years.csv", check.names = F)
creat <- read.csv("../Epigenetic Clocks WHO/Corr_Creatinine.csv")
creat = creat[,c(1,5)]
d1 = merge(d1, creat, by.x = "Sample_Name", by.y = "id", all.x = T)
names(d1)[138] <- "Creat"

females = d1[which(d1$sex %in% "F"),]
males = d1[which(d1$sex %in% "M"),] 

females_low <- females[which(females$Creat <= 62),] 
females_high <- females[which(females$Creat > 62),] 
females_na <- females[which(females$Creat %in% NA),] 

males_low <- males[which(males$Creat <= 80),] 
males_high <- males[which(males$Creat > 80),] 
males_na <- males[which(males$Creat %in% NA),] 

females_low$eGFR <- 141*((females_low$Creat/61.9)^-0.329)*(0.993^females_low$Age)*1.018
females_high$eGFR <- 141*((females_high$Creat/61.9)^-1.209)*(0.993^females_high$Age)*1.018
males_low$eGFR <- 141*((males_low$Creat/79.6)^-0.411)*(0.993^males_low$Age)
males_high$eGFR <-141*((males_high$Creat/79.6)^-1.209)*(0.993^males_high$Age)
females_na$eGFR <- NA
males_na$eGFR <- NA 
females_1 = rbind(females_low,females_high)
females_1 = rbind(females_1, females_na)
males_1 = rbind(males_low,males_high)
males_1 = rbind(males_1, males_na)
all = rbind(males_1, females_1)
d2<-all
d2$CKD = ifelse(d2$eGFR < 60, 1, 0)

write.csv(d2, "censoring_10years_CKD.csv",row.names=F)


## Censoring Date - 5 year cut-off ###

d1 <- read.csv("censoring_5years.csv", check.names = F)
creat <- read.csv("../Epigenetic Clocks WHO/Corr_Creatinine.csv")
creat = creat[,c(1,5)]
d1 = merge(d1, creat, by.x = "Sample_Name", by.y = "id", all.x = T)
names(d1)[138] <- "Creat"

females = d1[which(d1$sex %in% "F"),]
males = d1[which(d1$sex %in% "M"),] 

females_low <- females[which(females$Creat <= 62),] 
females_high <- females[which(females$Creat > 62),] 
females_na <- females[which(females$Creat %in% NA),] 

males_low <- males[which(males$Creat <= 80),] 
males_high <- males[which(males$Creat > 80),] 
males_na <- males[which(males$Creat %in% NA),] 

females_low$eGFR <- 141*((females_low$Creat/61.9)^-0.329)*(0.993^females_low$Age)*1.018
females_high$eGFR <- 141*((females_high$Creat/61.9)^-1.209)*(0.993^females_high$Age)*1.018
males_low$eGFR <- 141*((males_low$Creat/79.6)^-0.411)*(0.993^males_low$Age)
males_high$eGFR <-141*((males_high$Creat/79.6)^-1.209)*(0.993^males_high$Age)
females_na$eGFR <- NA
males_na$eGFR <- NA 
females_1 = rbind(females_low,females_high)
females_1 = rbind(females_1, females_na)
males_1 = rbind(males_low,males_high)
males_1 = rbind(males_1, males_na)
all = rbind(males_1, females_1)
d2<-all
d2$CKD = ifelse(d2$eGFR < 60, 1, 0)

write.csv(d2, "censoring_5years_CKD.csv",row.names=F)
