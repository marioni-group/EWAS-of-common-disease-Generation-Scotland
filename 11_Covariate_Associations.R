### Covariate prediction of incident disease ###

## Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/")

## Load requisite libraries 
library(survival)
library(reshape) 
library(ggplot2)
library(ggpubr)

## Read in covariate information 
cov=read.csv("covariates.csv")
## Only include those whose units/week of alcohol consumed was typical 

## Extract column names for analyses 
cols=c("Bcell", "CD4T", "CD8T", "Gran", "NK", "Mono", "units", "smokingScore", "rank", "years", "bmi")
names(cov)[20:ncol(cov)]=paste0("PC", 1:20)
cols2=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")

## Extract prevalence phenotypes
prev.phen=list.files("Prevalent_Phenotypes/", ".phen")
prev.phen=prev.phen[which(!prev.phen %in% c("asthma.phen","parkinsons_family.phen"))]
prev.phen.sex=prev.phen[which(prev.phen %in% c("breast_cancer.phen","prostate_cancer.phen"))]
prev.phen=prev.phen[which(!prev.phen %in% c("breast_cancer.phen","prostate_cancer.phen"))]

## Prepare lists to store output 
out.uni.total=list()
out.multi.total=list()

######### NON SEX-SPECIFIC ANALYSES #####

## start loop
for(i in prev.phen){ 
# read in file
tmp=read.table(paste0("Prevalent_Phenotypes/",i),header=T)  
# merge in with covariate information
tmp1=merge(tmp,cov,by.x="FID",by.y="Sample_Sentrix_ID",all.x=T)
# Create lists 
out.uni=list()
out.multi=list()
# loop through the covariates 
## Cells
for(j in c(cols,cols2)){ 
## Run glm models
out.uni[[j]] = summary(glm(phen ~ scale(tmp1[,j]), data = tmp1, family = "binomial"))$coefficients[2,]
out.multi[[j]] = summary(glm(phen ~ scale(tmp1[,j]) + age + sex, data= tmp1, family = "binomial"))$coefficients[2,]
}
## Tidy up dataframes 
out.uni=as.data.frame(do.call("rbind", out.uni))
out.multi=as.data.frame(do.call("rbind", out.multi))
out.uni$trait <- row.names(out.uni)
out.multi$trait <- row.names(out.multi)
out.uni$disease <- gsub(".phen.*", "", i)
out.multi$disease <- gsub(".phen.*", "", i)
out.uni.total[[i]] <- out.uni
out.multi.total[[i]] <- out.multi 
}

## Tidy up files so far 
uni.total=as.data.frame(do.call("rbind", out.uni.total))
multi.total=as.data.frame(do.call("rbind", out.multi.total))


######### SEX-SPECIFIC ANALYSES #####

## Prepare lists to store output 
out.uni.total=list()
out.multi.total=list()

## start loop
for(i in prev.phen.sex){ 
  # read in file
  tmp=read.table(paste0("Prevalent_Phenotypes/",i),header=T)  
  # merge in with covariate information
  tmp1=merge(tmp,cov,by.x="FID",by.y="Sample_Sentrix_ID",all.x=T)
  # Create lists 
  out.uni=list()
  out.multi=list()
  # loop through the covariates 
  ## Cells
  for(j in c(cols,cols2)){ 
    ## Run glm models
    out.uni[[j]] = summary(glm(phen ~ scale(tmp1[,j]), data = tmp1, family = "binomial"))$coefficients[2,]
    out.multi[[j]] = summary(glm(phen ~ scale(tmp1[,j]) + age, data= tmp1, family = "binomial"))$coefficients[2,]
  }
  ## Tidy up dataframes 
  out.uni=as.data.frame(do.call("rbind", out.uni))
  out.multi=as.data.frame(do.call("rbind", out.multi))
  out.uni$trait <- row.names(out.uni)
  out.multi$trait <- row.names(out.multi)
  out.uni$disease <- gsub(".phen.*", "", i)
  out.multi$disease <- gsub(".phen.*", "", i)
  out.uni.total[[i]] <- out.uni
  out.multi.total[[i]] <- out.multi 
}

## Tidy up files so far 
uni.total.sex=as.data.frame(do.call("rbind", out.uni.total))
multi.total.sex=as.data.frame(do.call("rbind", out.multi.total))

## Final files 
uni=rbind(uni.total, uni.total.sex)
multi=rbind(multi.total, multi.total.sex)

write.csv(uni, "Cleaned_Results/prevalent_univariate_covariate_associations.csv", row.names = F)
write.csv(multi, "Cleaned_Results/prevalent_multivariate_covariate_associations.csv", row.names = F)


## Extract incident phenotypes
inc.phen=list.files("Incident_Phenotypes/", ".phen")
inc.phen=inc.phen[which(!inc.phen %in% c("parkinsons_age65.phen","cervical_cancer.phen", "long_covid.phen","covid_hospitalisation.phen"))]
inc.phen.sex=inc.phen[which(inc.phen %in% c("breast_cancer.phen","prostate_cancer.phen","ovarian_cancer.phen"))]
inc.phen=inc.phen[which(!inc.phen %in% c("breast_cancer.phen","prostate_cancer.phen","ovarian_cancer.phen"))]


## Prepare lists to store output 
out.uni.total=list()
out.multi.total=list()

######### NON SEX-SPECIFIC ANALYSES #####

## start loop
for(i in inc.phen){ 
  # read in file
  tmp=read.table(paste0("Incident_Phenotypes/",i),header=T)  
  # merge in with covariate information
  tmp1=merge(tmp,cov,by="Sample_Sentrix_ID",all.x=T)
  # Create lists 
  out.uni=list()
  out.multi=list()
  # loop through the covariates 
  ## Cells
  for(j in c(cols,cols2)){ 
    ## Run glm models
    out.uni[[j]] = summary(coxph(Surv(tte,Event) ~ scale(tmp1[,j]), data = tmp1))$coefficients[1,]
    out.multi[[j]] = summary(coxph(Surv(tte,Event) ~ scale(tmp1[,j]) + age + sex, data = tmp1))$coefficients[1,]
  }
  ## Tidy up dataframes 
  out.uni=as.data.frame(do.call("rbind", out.uni))
  out.multi=as.data.frame(do.call("rbind", out.multi))
  out.uni$trait <- row.names(out.uni)
  out.multi$trait <- row.names(out.multi)
  out.uni$disease <- gsub(".phen.*", "", i)
  out.multi$disease <- gsub(".phen.*", "", i)
  out.uni.total[[i]] <- out.uni
  out.multi.total[[i]] <- out.multi 
}

## Tidy up files so far 
uni.total=as.data.frame(do.call("rbind", out.uni.total))
multi.total=as.data.frame(do.call("rbind", out.multi.total))


######### SEX-SPECIFIC ANALYSES #####

## Prepare lists to store output 
out.uni.total=list()
out.multi.total=list()

## start loop
for(i in inc.phen.sex){ 
  # read in file
  tmp=read.table(paste0("Incident_Phenotypes/",i),header=T)  
  # merge in with covariate information
  tmp1=merge(tmp,cov,by="Sample_Sentrix_ID",all.x=T)
  # Create lists 
  out.uni=list()
  out.multi=list()
  # loop through the covariates 
  ## Cells
  for(j in c(cols,cols2)){ 
    ## Run glm models
    out.uni[[j]] = summary(coxph(Surv(tte,Event) ~ scale(tmp1[,j]), data = tmp1))$coefficients[1,]
    out.multi[[j]] = summary(coxph(Surv(tte,Event) ~ scale(tmp1[,j]) + age, data = tmp1))$coefficients[1,]
  }
  ## Tidy up dataframes 
  out.uni=as.data.frame(do.call("rbind", out.uni))
  out.multi=as.data.frame(do.call("rbind", out.multi))
  out.uni$trait <- row.names(out.uni)
  out.multi$trait <- row.names(out.multi)
  out.uni$disease <- gsub(".phen.*", "", i)
  out.multi$disease <- gsub(".phen.*", "", i)
  out.uni.total[[i]] <- out.uni
  out.multi.total[[i]] <- out.multi 
}

## Tidy up files so far 
uni.total.sex=as.data.frame(do.call("rbind", out.uni.total))
multi.total.sex=as.data.frame(do.call("rbind", out.multi.total))

## Final files 
uni.first=rbind(uni.total, uni.total.sex)
multi.first=rbind(multi.total, multi.total.sex)


## Extract incident phenotypes - binary outcome not tte 
inc.phen=list.files("Incident_Phenotypes/", ".phen")
inc.phen=inc.phen[which(inc.phen %in% c("long_covid.phen","covid_hospitalisation.phen"))]


## Prepare lists to store output 
out.uni.total=list()
out.multi.total=list()

######### NON SEX-SPECIFIC ANALYSES #####

## start loop
for(i in inc.phen){ 
  # read in file
  tmp=read.table(paste0("Incident_Phenotypes/",i),header=T)  
  # merge in with covariate information
  tmp1=merge(tmp,cov,by="Sample_Sentrix_ID",all.x=T)
  # Create lists 
  out.uni=list()
  out.multi=list()
  # loop through the covariates 
  ## Cells
  for(j in c(cols,cols2)){ 
    ## Run glm models
    out.uni[[j]] = summary(glm(Event ~ scale(tmp1[,j]), data = tmp1, family = "binomial"))$coefficients[2,]
    out.multi[[j]] = summary(glm(Event ~ scale(tmp1[,j]) + covidAge + sex.x, data= tmp1, family = "binomial"))$coefficients[2,]
  }
  ## Tidy up dataframes 
  out.uni=as.data.frame(do.call("rbind", out.uni))
  out.multi=as.data.frame(do.call("rbind", out.multi))
  out.uni$trait <- row.names(out.uni)
  out.multi$trait <- row.names(out.multi)
  out.uni$disease <- gsub(".phen.*", "", i)
  out.multi$disease <- gsub(".phen.*", "", i)
  out.uni.total[[i]] <- out.uni
  out.multi.total[[i]] <- out.multi 
}

## Tidy up files so far 
uni.total=as.data.frame(do.call("rbind", out.uni.total))
multi.total=as.data.frame(do.call("rbind", out.multi.total))

# Prep final files
uni.first[,2]=NULL
multi.first[,2]=NULL
names(uni.first)=names(uni.total)
names(multi.first)=names(multi.total)

## Final files 
uni=rbind(uni.first, uni.total)
multi=rbind(multi.first, multi.total)

write.csv(uni, "Cleaned_Results/incident_univariate_covariate_associations.csv", row.names = F)
write.csv(multi, "Cleaned_Results/incident_multivariate_covariate_associations.csv", row.names = F)



###### PLOTTING STAGE ######
prev.uni=read.csv("Cleaned_Results/prevalent_univariate_covariate_associations.csv")
prev.multi=read.csv("Cleaned_Results/prevalent_multivariate_covariate_associations.csv")
inc.uni=read.csv("Cleaned_Results/incident_univariate_covariate_associations.csv")
inc.multi=read.csv("Cleaned_Results/incident_multivariate_covariate_associations.csv")

## add significance stars 
prev.uni$stars <- cut(prev.uni[,4], breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
prev.multi$stars <- cut(prev.multi[,4], breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
inc.uni$stars <- cut(inc.uni[,4], breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
inc.multi$stars <- cut(inc.multi[,4], breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

# Set full disease names 
prev.uni=prev.uni[order(prev.uni$disease),]
prev.multi=prev.multi[order(prev.multi$disease),]
inc.uni=inc.uni[order(inc.uni$disease),]
inc.multi=inc.multi[order(inc.multi$disease),]

prev.dis=c("Alzheimer's Disease", "Colorectal Cancer", "Breast Cancer", "Chronic Kidney Disease", "COPD", "Diabetes (Type 2)", "Ischemic Heart Disease", "Lung Cancer", "Osteoarthritis", "Chronic Neck/Back Pain", "Parkinson's Disease", "Prostate Cancer", "Rheumatoid Arthritis", "Stroke")
inc.dis=c("Alzheimer's Disease", "Colorectal Cancer", "Breast Cancer", "Chronic Kidney Disease", "COPD", "COVID - Hospitalisation", "Diabetes (Type 2)", "Ischemic Heart Disease", "Inflammatory Bowel Disease", "Liver Cirrhosis", "Long COVID", "Lung Cancer", "Osteoarthritis", "Ovarian Cancer", "Chronic Neck/Back Pain", "Parkinson's Disease", "Prostate Cancer", "Rheumatoid Arthritis", "Stroke")
prev.dis1=rep(prev.dis,each=31)
inc.dis1=rep(inc.dis,each=31)
prev.uni$disease=prev.dis1
prev.multi$disease=prev.dis1
inc.uni$disease=inc.dis1
inc.multi$disease=inc.dis1


## Reorder phenotypes 
prev.uni=prev.uni[order(prev.uni$disease),]
prev.multi=prev.multi[order(prev.multi$disease),]
inc.uni=inc.uni[order(inc.uni$disease),]
inc.multi=inc.multi[order(inc.multi$disease),]

prev.uni=prev.uni[which(prev.uni$trait %in% c("PC20", "PC19", "PC18", "PC17", "PC16", "PC15", "PC14", "PC13", "PC12", "PC11", "PC10","PC9","PC8","PC7","PC6","PC5","PC4","PC3","PC2","PC1","units", "bmi", "years", "rank", "smokingScore", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK")),]
prev.multi=prev.multi[which(prev.multi$trait %in% c("PC20", "PC19", "PC18", "PC17", "PC16", "PC15", "PC14", "PC13", "PC12", "PC11", "PC10","PC9","PC8","PC7","PC6","PC5","PC4","PC3","PC2","PC1","units", "bmi", "years", "rank", "smokingScore", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK")),]
inc.uni=inc.uni[which(inc.uni$trait %in% c("PC20", "PC19", "PC18", "PC17", "PC16", "PC15", "PC14", "PC13", "PC12", "PC11", "PC10","PC9","PC8","PC7","PC6","PC5","PC4","PC3","PC2","PC1","units", "bmi", "years", "rank", "smokingScore", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK")),]
inc.multi=inc.multi[which(inc.multi$trait %in% c("PC20", "PC19", "PC18", "PC17", "PC16", "PC15", "PC14", "PC13", "PC12", "PC11", "PC10","PC9","PC8","PC7","PC6","PC5","PC4","PC3","PC2","PC1","units", "bmi", "years", "rank", "smokingScore", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK")),]
 
prev.uni$trait <- factor(prev.uni$trait, levels = c("PC20", "PC19", "PC18", "PC17", "PC16", "PC15", "PC14", "PC13", "PC12", "PC11", "PC10","PC9","PC8","PC7","PC6","PC5","PC4","PC3","PC2","PC1","units", "bmi", "years", "rank", "smokingScore", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK"))
prev.multi$trait <- factor(prev.multi$trait, levels = c("PC20", "PC19", "PC18", "PC17", "PC16", "PC15", "PC14", "PC13", "PC12", "PC11", "PC10","PC9","PC8","PC7","PC6","PC5","PC4","PC3","PC2","PC1","units", "bmi", "years", "rank", "smokingScore", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK"))
inc.uni$trait <- factor(inc.uni$trait, levels = c("PC20", "PC19", "PC18", "PC17", "PC16", "PC15", "PC14", "PC13", "PC12", "PC11", "PC10","PC9","PC8","PC7","PC6","PC5","PC4","PC3","PC2","PC1","units", "bmi", "years", "rank", "smokingScore", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK"))
inc.multi$trait <- factor(inc.multi$trait, levels = c("PC20", "PC19", "PC18", "PC17", "PC16", "PC15", "PC14", "PC13", "PC12", "PC11", "PC10","PC9","PC8","PC7","PC6","PC5","PC4","PC3","PC2","PC1","units", "bmi", "years", "rank", "smokingScore", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK"))

##### PREVALENCE PLOTS #########

## PREVALENT UNIVARIATE PLOT 
p <- ggplot(aes(x=disease, y=trait, fill=z.value), data=prev.uni)
prev.uni.plot <- p + geom_tile() + scale_fill_gradient2(high="#D7191C", mid="white", low="#2C7BB6") + 
  geom_text(aes(label=stars), color="black", size=5) + 
  labs(y=NULL, x=NULL, fill="Z-score") + scale_y_discrete(labels=c("PC20", "PC19", "PC18", "PC17", "PC16", "PC15", "PC14", "PC13", "PC12", "PC11", "PC10","PC9","PC8","PC7","PC6","PC5","PC4","PC3","PC2","PC1","Alcohol", "Body Mass Index", "Education", "SIMD", "SmokingScore", "B cells", "CD4+T", "CD8+T","Granulocytes","Monocytes", "NK Cells")) +
  theme_bw() + ggtitle("Prevalent - Univariate Associations")+ theme(axis.text.x=element_text(angle = -60, hjust = 0),plot.title = element_text(hjust = 0.5)) 

## PREVALENT MULTIVARIATE PLOT 
p <- ggplot(aes(x=disease, y=trait, fill=z.value), data=prev.multi)
prev.multi.plot <- p + geom_tile() + scale_fill_gradient2(high="#D7191C", mid="white", low="#2C7BB6") + 
  geom_text(aes(label=stars), color="black", size=5) + 
  labs(y=NULL, x=NULL, fill="Z-score")  + scale_y_discrete(labels=c("PC20", "PC19", "PC18", "PC17", "PC16", "PC15", "PC14", "PC13", "PC12", "PC11", "PC10","PC9","PC8","PC7","PC6","PC5","PC4","PC3","PC2","PC1","Alcohol", "Body Mass Index", "Education", "SIMD", "SmokingScore", "B cells", "CD4+T", "CD8+T","Granulocytes","Monocytes", "NK Cells")) +
  theme_bw() + ggtitle("Prevalent - Multivariate Associations")+ theme(axis.text.x=element_text(angle = -60, hjust = 0),plot.title = element_text(hjust = 0.5))

pdf("Cleaned_Results/prevalent_covariates.pdf", width=14, height = 11, onefile=FALSE)
print(ggarrange(prev.uni.plot, prev.multi.plot, labels = c("A", "B"),align='h', common.legend=T, legend="right"))
dev.off()



##### INCIDENCE PLOTS #########

## INCIDENT UNIVARIATE PLOT 
## PREVALENT UNIVARIATE PLOT 
p <- ggplot(aes(x=disease, y=trait, fill=z.value), data=inc.uni)
inc.uni.plot <- p + geom_tile() + scale_fill_gradient2(high="#D7191C", mid="white", low="#2C7BB6") + 
  geom_text(aes(label=stars), color="black", size=5) + 
  labs(y=NULL, x=NULL, fill="Z-score") + scale_y_discrete(labels=c("PC20", "PC19", "PC18", "PC17", "PC16", "PC15", "PC14", "PC13", "PC12", "PC11", "PC10","PC9","PC8","PC7","PC6","PC5","PC4","PC3","PC2","PC1","Alcohol", "Body Mass Index", "Education", "SIMD", "SmokingScore", "B cells", "CD4+T", "CD8+T","Granulocytes","Monocytes", "NK Cells")) +
  theme_bw() + ggtitle("Incident - Univariate Associations")+ theme(axis.text.x=element_text(angle = -60, hjust = 0),plot.title = element_text(hjust = 0.5)) 

## PREVALENT MULTIVARIATE PLOT 
p <- ggplot(aes(x=disease, y=trait, fill=z.value), data=inc.multi)
inc.multi.plot <- p + geom_tile() + scale_fill_gradient2(high="#D7191C", mid="white", low="#2C7BB6") + 
  geom_text(aes(label=stars), color="black", size=5) + 
  labs(y=NULL, x=NULL, fill="Z-score")  + scale_y_discrete(labels=c("PC20", "PC19", "PC18", "PC17", "PC16", "PC15", "PC14", "PC13", "PC12", "PC11", "PC10","PC9","PC8","PC7","PC6","PC5","PC4","PC3","PC2","PC1","Alcohol", "Body Mass Index", "Education", "SIMD", "SmokingScore", "B cells", "CD4+T", "CD8+T","Granulocytes","Monocytes", "NK Cells")) +
  theme_bw() + ggtitle("Incident - Multivariate Associations")+ theme(axis.text.x=element_text(angle = -60, hjust = 0),plot.title = element_text(hjust = 0.5))

pdf("Cleaned_Results/incident_covariates.pdf", width=14, height = 11, onefile=FALSE)
print(ggarrange(inc.uni.plot, inc.multi.plot, labels = c("A", "B"),align='h', common.legend=T, legend="right"))
dev.off()

