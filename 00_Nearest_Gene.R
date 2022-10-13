# Load requisite libraries 
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)

# Prepare data 
# Get gene ids and their symbols 
annos=as.data.frame(transcripts(Homo.sapiens, columns=c("ENTREZID", "SYMBOL")))
genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
a=as.data.frame(genes)

####### INCIDENCE RESULTS ########

# Read in CpGs without annotations 
inc=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/incident_full_withreplication.csv")
# Set up dataframe to store outputs 
res=as.data.frame(matrix(nrow=nrow(inc),ncol=2))
names(res)=c("Overlap","new_gene")

# Loop stage
for(i in 1:nrow(inc)){ 
# Extract key info
chr=paste0("chr",inc[i,"Chr"])
bp=inc[i,"bp"]
# Query the chromosome and position
tmp=as.data.frame(genes[nearest(GRanges(chr, IRanges(bp,bp)),genes),])
# Get gene symbol
ans=unique(annos[which(annos$ENTREZID %in% tmp[,"gene_id"]),"SYMBOL"])[[1]]
# Store the output 
res[i,1]=inc[i,"Overlap"]
res[i,2]=ans
# Print to denote completion
print(i)
print(length(unique(annos[which(annos$ENTREZID %in% tmp[,"gene_id"]),"SYMBOL"])))
}

# Merge in with initial dataset and confirm that the annotation has worked 
res.inc=merge(inc,res,by="Overlap")

write.csv(res.inc, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/incident_full_replication_annotations.csv",row.names=F)



####### PREVALENCE RESULTS ########

# Read in CpGs without annotations 
inc=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/prevalent_full_withreplication.csv")

# Set up dataframe to store outputs 
res=as.data.frame(matrix(nrow=nrow(inc),ncol=2))
names(res)=c("Overlap","new_gene")

# Loop stage
for(i in 1:nrow(inc)){ 
  # Extract key info
  chr=paste0("chr",inc[i,"Chr"])
  bp=inc[i,"bp"]
  # Query the chromosome and position
  tmp=as.data.frame(genes[nearest(GRanges(chr, IRanges(bp,bp)),genes),])
  # Get gene symbol
  ans=unique(annos[which(annos$ENTREZID %in% tmp[,"gene_id"]),"SYMBOL"])[[1]]
  # Store the output 
  res[i,1]=inc[i,"Overlap"]
  res[i,2]=ans
  # Print to denote completion
  print(i)
  print(length(unique(annos[which(annos$ENTREZID %in% tmp[,"gene_id"]),"SYMBOL"])))
}

# Merge in with initial dataset and confirm that the annotation has worked 
res.inc=merge(inc,res,by="Overlap")
# res.inc[which(is.na(res.inc$new_gene)),"new_gene"]="ERICH1-AS1"

write.csv(res.inc, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/prevalent_full_replication_annotations.csv",row.names=F)