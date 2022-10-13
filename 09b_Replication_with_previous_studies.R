## Read in our results 
vars=c("prevalent_full","prevalent_wbcs","incident_full","incident_wbcs")
for(j in vars){ 
query=read.csv(paste0("../../../Results/", j, ".csv"))
query$Replicated=0
## Determine degree of replication within the existing studies
## extract files
files1=list.files(".", ".csv")
lit=gsub(".csv", "", files1)

# Identify which traits can be checked 
query1=query[which(query$trait %in% lit),]
query2=query[which(!query$trait %in% lit),]

## create output df 
out=as.data.frame(matrix(nrow=length(unique(query1$trait)),ncol=9))
names(out)=c("trait", "no_studies", "no_unique_cpgs_our_study", "no_unique_cpgs_literature", "no_unique_genes_our_study", "no_unique_genes_literature", "no_replicated_cpgs","no_replicated_genes", "gene_list")

# define list to store outputs 
list1=list()
# loop through files 
for(i in unique(query1$trait)){ 
  # get index of trait in queue 
  ind=which(unique(query1$trait) %in% i)
  ## identify existing studies on same trait 
  tmp=read.csv(paste0(i, ".csv"))
  # Subset our results to this trait 
  tmp1=query1[query1$trait %in% i,]
  ## extract trait name
  out[ind,1] <- as.character(i)
  ## extract no. of studies 
  out[ind,2] <- length(unique(tmp$author))
  ## extract no. of unique cpgs in our study
  out[ind,3] <- length(unique(tmp1$Probe))
  ## extract no. of unique cpgs in existing studies
  out[ind,4] <- length(unique(tmp$CpG))
  ## extract no. of unique genes in our study
  a=unique(tmp1$Gene)
  a=a[!is.na(a)]
  a=gsub(";.*","",a)
  a=unique(a)
  out[ind,5] <- length(a)
  ## extract no. of unique genes in existing studies 
  a1=unique(tmp$Gene)
  a1=a1[!is.na(a1)]
  a1=unique(a1)
  out[ind,6] <- length(a1)
  ## extract no. of replicated cpgs 
  out[ind,7] <- length(which(unique(tmp1$Probe) %in% unique(tmp$CpG)))
  ## extract proportion of replicated cpgs 
  out[ind,8] <- length(which(a %in% a1))
  out[ind,9] <- paste(a[which(a %in% a1)],collapse=",")
  # identify which CpGs are replicated 
  tmp1[which(tmp1$Probe %in% tmp$CpG),"Replicated"]=1
  list1[[i]]=tmp1
}

# combine replication query results and add in with traits that could not be tested 
list2=do.call("rbind", list1)
tmp2=rbind(list2,query2)
# save out files 
write.csv(tmp2, paste0("../../../Results/",j, "_withreplication.csv"),row.names=F)
write.csv(out, paste0("../../../Results/replication_", j, ".csv"),row.names=F)
print(j)
print(paste0(length(unique(query1$trait)),"/",length(unique(query$trait))))
}