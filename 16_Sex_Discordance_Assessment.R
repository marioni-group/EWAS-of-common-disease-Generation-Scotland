## Set working directory - Incident Overlap first 
setwd("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident_Sex/") 

# set up groups for analyses
split=c("males","females")
model=c("wbcs","PCs")

# Set up output df
output=as.data.frame(matrix(nrow=length(model),ncol=7))
names(output)=c("model","nCpG-males","nCpG-females","percent_agreed_males","percent_agreed_females","percent_agreed_males_1e5","percent_agreed_females_1e5")
for(j in 1:length(model)){ 
# males only
files=list.files(paste0(split[1],"/",model[j]),".linear")
files=files[which(!files%in%c("parkinsons_age65.linear","breast_cancer.linear","ovarian_cancer.linear", "prostate_cancer.linear", "cervical_cancer.linear"))]
files1=paste0(paste0(split[1],"/",model[j],"/"),files)
tmp=lapply(files1,fread)
names(tmp)=gsub(".linear", "", files)
tmp = Map(cbind, tmp,"trait"=names(tmp))
tmp1 <- as.data.frame(do.call(rbind, lapply(tmp, function(x)x[x$p < 1.9e-9,])))
tmp1.full <- as.data.frame(do.call(rbind, tmp))
tmp1.less <- as.data.frame(do.call(rbind, lapply(tmp, function(x)x[x$p < 0.05,])))

# females only
files2=list.files(paste0(split[2],"/",model[j]),".linear")
files2=files2[which(!files2%in%c("parkinsons_age65.linear","breast_cancer.linear","ovarian_cancer.linear","prostate_cancer.linear", "cervical_cancer.linear"))]
files3=paste0(paste0(split[2],"/",model[j],"/"),files2)
tmp2=lapply(files3,fread)
names(tmp2)=gsub(".linear", "", files2)
tmp2 = Map(cbind, tmp2,"trait"=names(tmp2))
tmp3 <- as.data.frame(do.call(rbind, lapply(tmp2, function(x)x[x$p < 1.9e-9,])))

# Determine overlap
tmp1$check=paste(tmp1$Probe,tmp1$trait,sep="_")
tmp1.less$check=paste(tmp1.less$Probe,tmp1.less$trait,sep="_")
tmp3$check=paste(tmp3$Probe,tmp3$trait,sep="_")
tmp3.less$check=paste(tmp3.less$Probe,tmp3.less$trait,sep="_")

# Store outputs
output[j,1]=as.character(model[[j]])
output[j,2]=nrow(tmp1)
output[j,3]=nrow(tmp3)
output[j,4]=signif(length(which(tmp1$check %in% tmp3$check))/nrow(tmp1),2)
output[j,5]=signif(length(which(tmp3$check %in% tmp1$check))/nrow(tmp3),2)
output[j,6]=signif(length(which(tmp1$check %in% tmp3.less$check))/nrow(tmp1),2)
output[j,7]=signif(length(which(tmp3$check %in% tmp1.less$check))/nrow(tmp3),2)

## Find discordant CpGs
tmp.m=tmp1[which(!tmp1$check%in%tmp3$check),]
tmp.f=tmp3[which(!tmp3$check%in%tmp1$check),]
tmp.m=tmp.m[order(tmp.m$p),]
tmp.f=tmp.f[order(tmp.f$p),]
tmp.m$group="males"
tmp.f$group="females"
tmp.disc=rbind(tmp.m,tmp.f)
## Save discordant outputs
write.csv(tmp.disc,paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Sex_Disagreements/Incident_",as.character(model[[j]]), ".csv"),row.names=F)


tmp=lapply(tmp,as.data.frame)
tmp2=lapply(tmp2,as.data.frame)
out.cor=as.data.frame(matrix(nrow=length(tmp),ncol=2))
for(i in 1:length(tmp)){
out.cor[i,1]=names(tmp)[i]
out.cor[i,2]=cor.test((tmp[[i]][,6]/tmp[[i]][,7]),(tmp2[[i]][,6]/tmp2[[i]][,7]))$estimate[[1]]
}
# # Plot agreement of beta values 
# tmp1.betas=tmp1[match(intersect(tmp1$check,tmp3$check),tmp1$check),]
# tmp3.betas=tmp3[match(intersect(tmp3$check,tmp1$check),tmp3$check),]
# Print to denote completion
print(j)
} 

# Write output 
write.csv(output, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident_Sex_overlap",row.names=F)


## Set working directory - Prevalent Overlap  
setwd("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent_Sex/") 

# set up groups for analyses
split=c("males","females")
model=c("wbcs","PCs")

# Set up output df
output=as.data.frame(matrix(nrow=length(model),ncol=7))
names(output)=c("model","nCpG-males","nCpG-females","percent_agreed_males","percent_agreed_females","percent_agreed_males_1e5","percent_agreed_females_1e5")
for(j in 1:length(model)){ 
  # males only
  files=list.files(paste0(split[1],"/",model[j]),".linear")
  files=files[which(!files%in%c("asthma.linear","breast_cancer.linear","parkinsons_family.linear","osteo_arthritis.linear","ovarian_cancer.linear", "prostate_cancer.linear", "cervical_cancer.linear"))]
  files1=paste0(paste0(split[1],"/",model[j],"/"),files)
  tmp=lapply(files1,fread)
  names(tmp)=gsub(".linear", "", files)
  tmp = Map(cbind, tmp,"trait"=names(tmp))
  tmp1 <- as.data.frame(do.call(rbind, lapply(tmp, function(x)x[x$p < 3.6e-8,])))
  tmp1.full <- as.data.frame(do.call(rbind, tmp))
  tmp1.less <- as.data.frame(do.call(rbind, lapply(tmp, function(x)x[x$p < 0.05,])))
  
  # females only
  files2=list.files(paste0(split[2],"/",model[j]),".linear")
  files2=files2[which(!files2%in%c("asthma.linear","parkinsons_family.linear","osteo_arthritis.linear","breast_cancer.linear","ovarian_cancer.linear","prostate_cancer.linear", "cervical_cancer.linear"))]
  files3=paste0(paste0(split[2],"/",model[j],"/"),files2)
  tmp2=lapply(files3,fread)
  names(tmp2)=gsub(".linear", "", files2)
  tmp2 = Map(cbind, tmp2,"trait"=names(tmp2))
  tmp3 <- as.data.frame(do.call(rbind, lapply(tmp2, function(x)x[x$p < 3.6e-8,])))
  tmp3.less <- as.data.frame(do.call(rbind, lapply(tmp2, function(x)x[x$p < 0.05,])))
  
  # Determine overlap
  tmp1$check=paste(tmp1$Probe,tmp1$trait,sep="_")
  tmp1.less$check=paste(tmp1.less$Probe,tmp1.less$trait,sep="_")
  tmp3$check=paste(tmp3$Probe,tmp3$trait,sep="_")
  tmp3.less$check=paste(tmp3.less$Probe,tmp3.less$trait,sep="_")
  
  # Store outputs
  output[j,1]=as.character(model[[j]])
  output[j,2]=nrow(tmp1)
  output[j,3]=nrow(tmp3)
  output[j,4]=signif(length(which(tmp1$check %in% tmp3$check))/nrow(tmp1),2)
  output[j,5]=signif(length(which(tmp3$check %in% tmp1$check))/nrow(tmp3),2)
  output[j,6]=signif(length(which(tmp1$check %in% tmp3.less$check))/nrow(tmp1),2)
  output[j,7]=signif(length(which(tmp3$check %in% tmp1.less$check))/nrow(tmp3),2)
  
  ## Find discordant CpGs
  tmp.m=tmp1[which(!tmp1$check%in%tmp3$check),]
  tmp.f=tmp3[which(!tmp3$check%in%tmp1$check),]
  tmp.m=tmp.m[order(tmp.m$p),]
  tmp.f=tmp.f[order(tmp.f$p),]
  tmp.m$group="males"
  tmp.f$group="females"
  tmp.disc=rbind(tmp.m,tmp.f)
  ## Save discordant outputs
  write.csv(tmp.disc,paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Sex_Disagreements/Prevalent_",as.character(model[[j]]), ".csv"),row.names=F)

  # # Plot agreement of beta values 
  # tmp1.betas=tmp1[match(intersect(tmp1$check,tmp3$check),tmp1$check),]
  # tmp3.betas=tmp3[match(intersect(tmp3$check,tmp1$check),tmp3$check),]
  # Print to denote completion
  print(j)
} 


tmp=lapply(tmp,as.data.frame)
tmp2=lapply(tmp2,as.data.frame)
out.cor=as.data.frame(matrix(nrow=length(tmp),ncol=2))
for(i in 1:length(tmp)){
  out.cor[i,1]=names(tmp)[i]
  out.cor[i,2]=cor.test((tmp[[i]][,6]/tmp[[i]][,7]),(tmp2[[i]][,6]/tmp2[[i]][,7]))$estimate[[1]]
}

# Write output 
write.csv(output, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent_Sex_overlap.csv",row.names=F)


## Read prevalent analyses 
out1=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Sex_Disagreements/correlations_prevalent.csv")
out1$Trait=c("Alzheimer's Disease", "Colorectal Cancer", "Chronic Kidney Disease", "COPD", "Diabetes (Type 2)", "Ischemic Heart Disease", "Lung Cancer", "Osteoarthritis", "Chronic Neck/Back Pain", "Parkinson's Disease", "Rhemumatoid Arthritis", "Stroke")
out1=out1[order(out1$Trait),]
names(out1)=c("Trait", "Basic","Full")
data <- gather(out1, condition, value, Basic:Full, factor_key=TRUE)
data$coef=signif(data$value,2)

p=ggplot(aes(x=Trait, y=condition, fill=value), data=data)
prev.uni.plot <- p + geom_tile() + scale_fill_gradient2(high="#D7191C", mid="white", low="#2C7BB6") + 
  geom_text(aes(label=coef), color="black", size=3.7) + 
  labs(y=NULL, x=NULL, fill="Correlation Coefficient between \n Z-scores from Male and Female-specific EWAS") + scale_y_discrete(labels=c("Basic Model","Fully-Adjusted Model"))+
  theme_bw() + ggtitle("Prevalence")+ theme(axis.text.y=element_text(size=10.5),axis.text.x=element_text(angle = -60, hjust = 0,size=10.5),plot.title = element_text(hjust = 0.5))  +theme(legend.position="none")



## Read incident analyses 


## Read prevalent analyses 
out2=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Sex_Disagreements/correlations_incident.csv")
out2$Trait=c("Alzheimer's Disease", "Colorectal Cancer", "Chronic Kidney Disease", "COPD", "Diabetes (Type 2)", "Ischemic Heart Disease", "Inflammatory Bowel Disease", "Liver Cirrhosis", "Lung Cancer", "Osteoarthritis", "Chronic Neck/Back Pain", "Parkinson's Disease", "Rhemumatoid Arthritis", "Stroke")
out2=out2[order(out2$Trait),]
names(out2)=c("Trait", "Basic","Full")
data1 <- gather(out2, condition, value, Basic:Full, factor_key=TRUE)
data1$coef=signif(data1$value,2)

p1=ggplot(aes(x=Trait, y=condition, fill=value), data=data1)
inc.uni.plot <- p1 + geom_tile() + scale_fill_gradient2(high="#D7191C", mid="white", low="#2C7BB6") + 
  geom_text(aes(label=coef), color="black", size=3.7) + 
  labs(y=NULL, x=NULL, fill="Correlation Coefficient between \n Z-scores from Male and Female-specific EWAS") + scale_y_discrete(labels=c("Basic Model","Fully-Adjusted Model"))+
  theme_bw() + ggtitle("Incidence")+ theme(axis.text.y=element_text(size=10.5),axis.text.x=element_text(angle = -60, hjust = 0,size=10.5),plot.title = element_text(hjust = 0.5)) +theme(legend.position="none")


l1=ggplot(aes(x=Trait, y=condition, fill=value), data=data1)
graph_legend <- l1 + geom_tile() + scale_fill_gradient2(high="#D7191C", mid="white", low="#2C7BB6") + 
  geom_text(aes(label=coef), color="black", size=3.7) + 
  labs(y=NULL, x=NULL, fill="Correlation Coefficient between \nZ-scores from Male and Female-specific EWAS") + scale_y_discrete(labels=c("Basic Model","Fully-Adjusted Model"))+
  theme_bw() + ggtitle("Incidence")+ theme(axis.text.y=element_text(size=10.5),axis.text.x=element_text(angle = -60, hjust = 0,size=10.5),plot.title = element_text(hjust = 0.5)) 

library(cowplot)

legend <- get_legend(graph_legend)
pdf("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/sex_specific_correlations.pdf",width=14,height=8.5)
plot = plot_grid(plot_grid(prev.uni.plot, inc.uni.plot, nrow = 2), legend, ncol = 2, rel_widths = c(0.75,0.3))
print(plot)
dev.off()