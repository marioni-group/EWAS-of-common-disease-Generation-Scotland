## Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/")

## Load requisite libraries
library(data.table)

## List files 
## Prevalent
prev.wbcs=list.files("Prevalent/Outputs/wbcs/", ".linear")
prev.full=list.files("Prevalent/Outputs/PCs/", ".linear")
## Incident
inc.wbcs=list.files("Incident/Outputs/wbcs/", ".linear")
inc.full=list.files("Incident/Outputs/PCs/", ".linear")
## Incident - 5 year
inc.5.wbcs=list.files("Incident_5years/Outputs/wbcs/", ".linear")
inc.5.full=list.files("Incident_5years/Outputs/PCs/", ".linear")
## Incident - 10 year
inc.10.wbcs=list.files("Incident_10years/Outputs/PCs/", ".linear")
inc.10.full=list.files("Incident_10years/Outputs/full/", ".linear")

## Make lists of all significant CpGs (P<3.6e-8 - can filter to Bonferroni threshold later)

## Set up lists to store outputs 
lists=list(prev.wbcs,prev.full,inc.wbcs,inc.full,inc.5.wbcs,inc.5.full,inc.10.wbcs,inc.10.full)
list.paths=c("Prevalent/Outputs/wbcs/","Prevalent/Outputs/PCs/", "Incident/Outputs/wbcs/" ,"Incident/Outputs/PCs/", "Incident_5years/Outputs/wbcs/" ,"Incident_5years/Outputs/PCs/", "Incident_10years/Outputs/wbcs/" ,"Incident_10years/Outputs/PCs/")
out.names=c("prevalent_wbcs_08102022", "prevalent_full_08102022", "incident_wbcs_08102022","incident_full_08102022", "incident_5years_wbcs_08102022","incident_5years_full_08102022", "incident_10years_wbcs_08102022","incident_10years_full_08102022")

## start loop
for(j in 1:length(lists)){
## Create final list to store CpGs 
out<-list()
## Loop through results files within individual folders
for(i in lists[[j]]){
## Read in file
tmp=as.data.frame(fread(paste0(list.paths[j],i)))  
## Store trait name
tmp$trait <- gsub(".linear.*", "", i)
## Subset to signficant CpGs only 
tmp1=tmp[which(tmp$p < 3.6e-8),]
## If there are no CpGs, then skip 
if(nrow(tmp1)==0){ 
  NULL 
} ## If there are CpGs, then store them in list 
else { 
out[[i]] <- tmp1
}
## Print to denote completion
print(i)
}
## Tidy the resulting dataframe of CpGs for given folder 
out1 <- as.data.frame(do.call("rbind", out))
out1 <- out1[order(out1$trait, out1$p),]
## Write output 
write.csv(out1, paste0("Results/", out.names[j],"_output.csv"), row.names = F)
## Print to denote completion
print(j)
}


## Read in compiled results list and tidy 
setwd("Results/")
prev.wbcs=read.csv("prevalent_wbcs_08102022_output.csv")
prev.full=read.csv("prevalent_full_08102022_output.csv")
inc.wbcs=read.csv("incident_wbcs_08102022_output.csv")
inc.full=read.csv("incident_full_08102022_output.csv")
inc.5.wbcs=read.csv("incident_5years_wbcs_08102022_output.csv")
inc.5.full=read.csv("incident_5years_full_08102022_output.csv")
inc.10.wbcs=read.csv("incident_10years_wbcs_08102022_output.csv")
inc.10.full=read.csv("incident_10years_full_08102022_output.csv")

## Thresholds 
## Prevalent - 14 diseases at present
t1=3.6e-8/14
## Incident - 11 diseases at present 
t2=3.6e-8/19

## Read in 'good' probes 
keep=read.table("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt",header=F)

## Tidy up
tidy.prev <- function(df) {
  tmp1 <- df[which(df$p < t1),]
  tmp1 <- tmp1[which(tmp1$Probe %in% keep$V1),]
  tmp1 <- tmp1[which(!tmp1$trait %in% c("asthma","osteo_arthritis", "parkinsons_family", "parkinsons_age65", "cervical_cancer")),]
  return(tmp1)
}

tidy.inc <- function(df) {
  tmp1 <- df[which(df$p < t2),]
  tmp1 <- tmp1[which(tmp1$Probe %in% keep$V1),]
  tmp1 <- tmp1[which(!tmp1$trait %in% c("asthma","osteo_arthritis", "parkinsons_family", "parkinsons_age65", "cervical_cancer")),]
  return(tmp1)
}

prev.wbcs=tidy.prev(prev.wbcs)
prev.full=tidy.prev(prev.full)
inc.wbcs=tidy.inc(inc.wbcs)
inc.full=tidy.inc(inc.full)
inc.5.wbcs=tidy.inc(inc.5.wbcs)
inc.5.full=tidy.inc(inc.5.full)
inc.10.wbcs=tidy.inc(inc.10.wbcs)
inc.10.full=tidy.inc(inc.10.full)


## Count number of results 
nrow(prev.wbcs)
nrow(prev.full)
nrow(inc.wbcs)
nrow(inc.full)
nrow(inc.5.wbcs)
nrow(inc.5.full)
nrow(inc.10.wbcs)
nrow(inc.10.full)


## Compare overlap 
prev.wbcs$Overlap=paste(prev.wbcs$Probe, prev.wbcs$trait,sep="_")
prev.full$Overlap=paste(prev.full$Probe, prev.full$trait,sep="_")
inc.wbcs$Overlap=paste(inc.wbcs$Probe, inc.wbcs$trait,sep="_")
inc.full$Overlap=paste(inc.full$Probe, inc.full$trait,sep="_")
inc.5.wbcs$Overlap=paste(inc.5.wbcs$Probe, inc.5.wbcs$trait,sep="_")
inc.5.full$Overlap=paste(inc.5.full$Probe, inc.5.full$trait,sep="_")
inc.10.wbcs$Overlap=paste(inc.10.wbcs$Probe, inc.10.wbcs$trait,sep="_")
inc.10.full$Overlap=paste(inc.10.full$Probe, inc.10.full$trait,sep="_")


## compare basic models
prev.comp=prev.wbcs[which(prev.wbcs$trait %in% inc.wbcs$trait),]
inc.comp=inc.wbcs[which(inc.wbcs$trait %in% prev.wbcs$trait),]

length(which(prev.comp$Overlap %in% inc.comp$Overlap))/nrow(prev.comp)
length(which(inc.comp$Overlap %in% prev.comp$Overlap))/nrow(inc.comp)

## compare full models 
prev.comp1=prev.full[which(prev.full$trait %in% inc.full$trait),]
inc.comp1=inc.full[which(inc.full$trait %in% prev.full$trait),]

length(which(prev.comp1$Overlap %in% inc.comp1$Overlap))/nrow(prev.comp1)
length(which(inc.comp1$Overlap %in% prev.comp1$Overlap))/nrow(inc.comp1)


## compare agreement within prev and inc 
length(which(prev.full$Overlap %in% prev.wbcs$Overlap))/nrow(prev.full)
length(which(inc.full$Overlap %in% inc.wbcs$Overlap))/nrow(inc.full)

## Add extra columns 
prev.wbcs$Present_in_Full_Model=0
prev.wbcs$Present_in_Incidence_Model=0
prev.wbcs[which(prev.wbcs$Overlap %in% prev.full$Overlap),"Present_in_Full_Model"]=1
prev.wbcs[which(prev.wbcs$Overlap %in% inc.wbcs$Overlap),"Present_in_Incidence_Model"]=1

prev.full$Present_in_Basic_Model=0
prev.full$Present_in_Incidence_Model=0
prev.full[which(prev.full$Overlap %in% prev.wbcs$Overlap),"Present_in_Basic_Model"]=1
prev.full[which(prev.full$Overlap %in% inc.full$Overlap),"Present_in_Incidence_Model"]=1

inc.wbcs$Present_in_Full_Model=0
inc.wbcs$Present_in_Prevalence_Model=0
inc.wbcs$Present_at_5_Year_Follow_Up=0
inc.wbcs$Present_at_10_Year_Follow_Up=0
inc.wbcs[which(inc.wbcs$Overlap %in% inc.full$Overlap),"Present_in_Full_Model"]=1
inc.wbcs[which(inc.wbcs$Overlap %in% prev.wbcs$Overlap),"Present_in_Prevalence_Model"]=1
inc.wbcs[which(inc.wbcs$Overlap %in% inc.5.wbcs$Overlap),"Present_at_5_Year_Follow_Up"]=1
inc.wbcs[which(inc.wbcs$Overlap %in% inc.10.wbcs$Overlap),"Present_at_10_Year_Follow_Up"]=1
inc.5.wbcs$Present_in_Full_Follow_Up=0
inc.10.wbcs$Present_in_Full_Follow_Up=0
inc.5.wbcs[which(inc.5.wbcs$Overlap %in% inc.wbcs$Overlap),"Present_in_Full_Follow_Up"]=1
inc.10.wbcs[which(inc.10.wbcs$Overlap %in% inc.wbcs$Overlap),"Present_in_Full_Follow_Up"]=1


inc.full$Present_in_Basic_Model=0
inc.full$Present_in_Prevalence_Model=0
inc.full$Present_at_5_Year_Follow_Up=0
inc.full$Present_at_10_Year_Follow_Up=0
inc.full[which(inc.full$Overlap %in% inc.wbcs$Overlap),"Present_in_Basic_Model"]=1
inc.full[which(inc.full$Overlap %in% prev.full$Overlap),"Present_in_Prevalence_Model"]=1
inc.full[which(inc.full$Overlap %in% inc.5.full$Overlap),"Present_at_5_Year_Follow_Up"]=1
inc.full[which(inc.full$Overlap %in% inc.10.full$Overlap),"Present_at_10_Year_Follow_Up"]=1
inc.5.full$Present_in_Full_Follow_Up=0
inc.10.full$Present_in_Full_Follow_Up=0
inc.5.full[which(inc.5.full$Overlap %in% inc.full$Overlap),"Present_in_Full_Follow_Up"]=1
inc.10.full[which(inc.10.full$Overlap %in% inc.full$Overlap),"Present_in_Full_Follow_Up"]=1


## Write out cleaned files 
write.csv(prev.wbcs,"../Cleaned_Results/prevalent_wbcs.csv",row.names=F)
write.csv(prev.full,"../Cleaned_Results/prevalent_full.csv",row.names=F)
write.csv(inc.wbcs,"../Cleaned_Results/incident_wbcs.csv",row.names=F)
write.csv(inc.full,"../Cleaned_Results/incident_full.csv",row.names=F)
write.csv(inc.5.wbcs,"../Cleaned_Results/incident_5years_wbcs.csv",row.names=F)
write.csv(inc.5.full,"../Cleaned_Results/incident_5years_full.csv",row.names=F)
write.csv(inc.10.wbcs,"../Cleaned_Results/incident_10years_wbcs.csv",row.names=F)
write.csv(inc.10.full,"../Cleaned_Results/incident_10years_full.csv",row.names=F)
