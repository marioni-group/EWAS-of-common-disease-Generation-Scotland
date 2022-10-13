## Set working directory - Incident Overlap first 
setwd("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident_Sex/") 

# set up groups for analyses
split=c("males","females")
model=c("basic","wbcs","PCs")

# Set up output df
output=as.data.frame(matrix(nrow=length(model),ncol=7))
names(output)=c("model","nCpG-males","nCpG-females","percent_agreed_males","percent_agreed_females","percent_agreed_males_1e5","percent_agreed_females_1e5")
for(j in 1:length(model)){ 
# males only
files=list.files(paste0(split[1],"/",model[j]),".linear")
files=files[which(!files%in%c("parkinsons.linear","breast_cancer.linear","ovarian_cancer.linear", "prostate_cancer.linear", "cervical_cancer.linear"))]
files1=paste0(paste0(split[1],"/",model[j],"/"),files)
tmp=lapply(files1,fread)
names(tmp)=gsub(".linear", "", files)
tmp = Map(cbind, tmp,"trait"=names(tmp))
tmp1 <- as.data.frame(do.call(rbind, lapply(tmp, function(x)x[x$p < 3.6e-8,])))
tmp1.full <- as.data.frame(do.call(rbind, tmp))
tmp1.less <- as.data.frame(do.call(rbind, lapply(tmp, function(x)x[x$p < 1e-5,])))

# females only
files2=list.files(paste0(split[2],"/",model[j]),".linear")
files2=files2[which(!files2%in%c("parkinsons.linear","breast_cancer.linear","ovarian_cancer.linear","prostate_cancer.linear", "cervical_cancer.linear"))]
files3=paste0(paste0(split[2],"/",model[j],"/"),files2)
tmp2=lapply(files3,fread)
names(tmp2)=gsub(".linear", "", files2)
tmp2 = Map(cbind, tmp2,"trait"=names(tmp2))
tmp3 <- as.data.frame(do.call(rbind, lapply(tmp2, function(x)x[x$p < 3.6e-8,])))
tmp3.less <- as.data.frame(do.call(rbind, lapply(tmp2, function(x)x[x$p < 1e-5,])))

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
model=c("basic","wbcs","PCs")

# Set up output df
output=as.data.frame(matrix(nrow=length(model),ncol=7))
names(output)=c("model","nCpG-males","nCpG-females","percent_agreed_males","percent_agreed_females","percent_agreed_males_1e5","percent_agreed_females_1e5")
for(j in 1:length(model)){ 
  # males only
  files=list.files(paste0(split[1],"/",model[j]),".linear")
  files=files[which(!files%in%c("asthma.linear","breast_cancer.linear","ovarian_cancer.linear", "prostate_cancer.linear", "cervical_cancer.linear"))]
  files1=paste0(paste0(split[1],"/",model[j],"/"),files)
  tmp=lapply(files1,fread)
  names(tmp)=gsub(".linear", "", files)
  tmp = Map(cbind, tmp,"trait"=names(tmp))
  tmp1 <- as.data.frame(do.call(rbind, lapply(tmp, function(x)x[x$p < 3.6e-8,])))
  tmp1.full <- as.data.frame(do.call(rbind, tmp))
  tmp1.less <- as.data.frame(do.call(rbind, lapply(tmp, function(x)x[x$p < 1e-5,])))
  
  # females only
  files2=list.files(paste0(split[2],"/",model[j]),".linear")
  files2=files2[which(!files2%in%c("asthma.linear","breast_cancer.linear","ovarian_cancer.linear","prostate_cancer.linear", "cervical_cancer.linear"))]
  files3=paste0(paste0(split[2],"/",model[j],"/"),files2)
  tmp2=lapply(files3,fread)
  names(tmp2)=gsub(".linear", "", files2)
  tmp2 = Map(cbind, tmp2,"trait"=names(tmp2))
  tmp3 <- as.data.frame(do.call(rbind, lapply(tmp2, function(x)x[x$p < 3.6e-8,])))
  tmp3.less <- as.data.frame(do.call(rbind, lapply(tmp2, function(x)x[x$p < 1e-5,])))
  
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

# Write output 
write.csv(output, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident_Sex_overlap",row.names=F)

## Follow up testing 
inc.basic=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Sex_Disagreements/Incident_basic.csv")
inc.basic.m=inc.basic[inc.basic$group%in%"males",]
inc.basic.f=inc.basic[inc.basic$group%in%"females",]
gene.m=unique(inc.basic.m$Gene)
gene.f=unique(inc.basic.f$Gene)
gene.m.d=gene.m[!gene.m%in%gene.f]
gene.f.d=gene.f[!gene.f%in%gene.m]
