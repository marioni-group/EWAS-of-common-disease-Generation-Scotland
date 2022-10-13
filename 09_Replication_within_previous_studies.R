### Compile existing studies to test replication 

## set working directory 
setwd("C:/Users/rhillar2/Desktop/Rob/EWAS_Disease_GS/Literature/Final")

## set variables to extract 
vars=c("heart_disease", "stroke", "diabetes", "breast_cancer", "bowel_cancer", "prostate_cancer", "lung_cancer", "COPD", "pain", "rheumatoid_arthritis", "osteo_arthritis", "ibd", "parkinsons", "ad", "covid_severity", "chronic_kidney_disease", "long_covid", "liver_cirrhosis", "ovarian_cancer")

## extract all files in vector
files=list.files(".", ".csv")

## loop through each variable and extract relevant files
## then rbind-ing them together to make final file for that trait

for(i in vars){ 
## extract all files with variable name 
tmp.files=files[grep(i, files)]
## read in and rbind the files 
tmp.df=as.data.frame(do.call("rbind", lapply(tmp.files, function(x) read.csv(x,header=T))))
## print number of unique studies - sense check 
print(length(unique(tmp.df$author)))
## write out file if there are > 1 rows 
if(!nrow(tmp.df)==0){ 
write.csv(tmp.df, paste0("Summary/",i,".csv"),row.names = F)
} else { NULL }
} 

## Determine degree of replication within the existing studies
setwd("Summary/")
## extract files
files1=list.files(".", ".csv")

## create output df 
out=as.data.frame(matrix(nrow=length(files1),ncol=9))
names(out)=c("trait", "no_studies", "no_unique_cpgs", "no_replicated_cpgs", "prop_replicated","no_unique_genes","no_replicated_genes","prop_replicated_genes", "gene_list")

## loop through files 
for(i in 1:length(files1)){ 
## read in files
tmp=read.csv(files1[i])
if(length(unique(tmp$author))<2){ NULL } else {
## extract trait name
out[i,1] <- as.character(gsub(".csv.*", "", files1[i]))
## extract no. of studies 
out[i,2] <- length(unique(tmp$author))
## extract no. of unique cpgs 
cpgs=as.data.frame(table(tmp$CpG))
out[i,3] <- nrow(cpgs)
## extract no. of replicated cpgs 
out[i,4] <- length(which(cpgs$Freq>1))
## extract proportion of replicated cpgs 
out[i,5] <- signif(length(which(cpgs$Freq>1))/nrow(cpgs),2)
## extract no. of unique and replicated genes
genes=table(tmp$Gene, tmp$author)
genes=apply(genes, 1, function(x) x<-ifelse(x>1, 1, x))
genes=as.data.frame(colSums(genes))
genes$Gene = row.names(genes)
genes=genes[which(!genes$Gene %in% c("-","")),]
names(genes)[1]="Freq"
out[i,6] <- nrow(genes)
out[i,7] <- length(which(genes$Freq>1))
## extract proportion of replicated genes 
out[i,8] <- signif(length(which(genes$Freq>1))/nrow(genes),2)
out[i,9] <- paste(as.character(genes[genes$Freq>1,"Gene"]),collapse=",")
## subset to pe-7
## print to denote completion
print(files1[i])
}
}

out=out[which(!is.na(out$trait)),]
write.csv(out, "../../../Replication_between_previous_studies.csv",row.names=F)
