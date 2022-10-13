## Preparation of Time-To-Event variables 

## Set working directory
setwd("C:/Users/rhillar2/Desktop/Rob/EWAS_Disease_GS/Incidence/")

## Read in incident phenotypes 
Diabetes=read.csv("Diabetes_combined.csv")

## Read in censoring info 
cens=read.csv("../censoring_oct2020.csv")

######################################################
#### Diseases with prevalent data at baseline ########
######################################################

## Set up loop 
my.list = list(Diabetes)
my.names = c("diabetes_Y")

for(i in 1:length(my.list)){
  tmp <- my.list[[i]]
  ## Remove Prevalent cases 
  prev.tmp = cens[-which(cens[,my.names[[i]]]%in% 1),]
  tmp1 = tmp[which(tmp$id %in% prev.tmp$Sample_Name),]
  for(ind in 1:14){ 
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
  cox[which(cox$tte > ind),"Event"]<-0
  cox[which(cox$tte > ind),"tte"]<-ind
  cox1<-cox[,c("Sample_Name","Sample_Sentrix_ID","Event","tte")]
  cox1<-cox1[which(!is.na(cox1$tte)),]
  
  ## Extract variable name and write out file
  var.name=gsub("_Y", "", my.names[[i]])
  write.table(cox1,paste0("../CoxPH_Violations/", var.name, "/year_", ind, ".phen"),row.names=F, sep=' ')
  
  ## Print to denote completion
  print(var.name)
  
}
} 




