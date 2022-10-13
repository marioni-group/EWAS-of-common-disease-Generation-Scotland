
## Load requisite libraries 
library(lumi)
library(limma)

#################################################################
### STEP 1. IDENTIFY CpGs AND PREP PHENOTYPE FILES ##############
#################################################################

# Read in results file 
tmp=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Results/Results_Final/incident_full.csv")

# Remove traits with <20 individuals for analyses 
tmp=tmp[!tmp$trait %in% c("cervical_cancer"),]

# Read in Methylation Data and Sample Information 
dnam=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/mvals.rds")
samps=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")

# Subset CpGs to those with significant associations with disease 
dnam1=dnam[which(row.names(dnam)%in%tmp$Probe),]

# Subset individuals to those with genetic data 
gs=read.table("/Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/GS20K_chr10_HRC.r1-1_nomono_I4_cpra.fam") 
samps1=samps[samps$Sample_Name %in% gs$V1,]
dnam1=dnam1[,which(colnames(dnam1)%in%samps1$Sample_Sentrix_ID)]

# Convert to beta values 
dnam1=m2beta(dnam1)

# Next we need to regress out covariates as in GoDMC protocol 
  #-age, sex, batch, genetic PCs, predicted smoking, predicted cell counts 
cov=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/covariates.csv")
cov1=cov[,c("Sample_Sentrix_ID","age","sex","Batch","Gran", "NK","Bcell","CD4T", "CD8T", "smokingScore", "V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17","V18","V19","V20","V21", "V22")]
cov2=cov1[complete.cases(cov1),]
dnam1=dnam1[,which(colnames(dnam1)%in%cov2$Sample_Sentrix_ID)]
ids=cov2$Sample_Sentrix_ID
dnam1=dnam1[,match(ids,colnames(dnam1))]

 # Residualisation step 
design.resid <- model.matrix(~as.factor(sex) + age + as.factor(Batch) + smokingScore + NK + Gran + Bcell + CD4T + CD8T + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22, data=cov1)
fit.resid <- limma::lmFit(dnam1, design.resid)
gc()
dnam2 <- limma::residuals.MArrayLM(fit.resid, dnam1)
dnam2 <- dnam2[!is.infinite(rowSums(dnam2)), ]
rm(fit.resid)
gc()


# Transpose and merge in Sample_Name (GWAS ID) information so that we can create fastGWA phenotype files 
dnam3=as.data.frame(apply(dnam2,1,scale))
row.names(dnam3)=colnames(dnam2)
dnam3$Sample_Sentrix_ID=row.names(dnam3)
dnam3=merge(dnam3,samps1[,c("Sample_Sentrix_ID","Sample_Name")],by="Sample_Sentrix_ID")

# Make phenotypes for fastGWA 
for(i in 2:(ncol(dnam3)-1)){ 
  tmp <- dnam3[,c(ncol(dnam3),ncol(dnam3), i)]
  names(tmp)[1:2]=c("FID","IID")
  print(i)
  write.table(tmp, file = paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/mQTLs/Inputs/",names(dnam3)[i], ".phen"), row.names=F, col.names = F, sep=' ', quote = F)
}

## Create keep file for samps 
keep=tmp[,c("FID", "IID")]
write.table(keep, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/mQTLs/keep_samps.txt",row.names=F,col.names=F,quote=F,sep=" ")


####################################
### STEP 2. RUN GWAS  ##############
####################################

/home/robert/gcta_1.93.2beta/gcta64 --mbfile  /Cluster_Filespace/Marioni_Group/Daniel/Methyl_Protein_GWAS/gs_mbfile.txt --make-grm --sparse-cutoff 0.05 --keep /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/mQTLS/keep_samps.txt --threads 10 --out /Cluster_Filespace/Marioni_Group/Rob/Somalogic/methid_mqtl_grm 


cd /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/mQTLs/Inputs/
  for i in *.phen 

do

/home/robert/gcta_1.93.2beta/gcta64 --mbfile /Cluster_Filespace/Marioni_Group/Daniel/Methyl_Protein_GWAS/gs_mbfile.txt \
--grm-sparse /Cluster_Filespace/Marioni_Group/Rob/Somalogic/methid_mqtl_grm \
--fastGWA-mlm \
--maf 0.01 \
--pheno /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/mQTLs/Inputs/$i  \
--threads 10 \
--out /Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/mQTLs/Outputs/$i

done

#####################################
#### PROCESS GWAS RESULTS  ##########
#####################################
library(data.table)
files = list.files("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/mQTLs/Outputs/", ".fastGWA")

stats <- list()
for(i in files){ 
  tmp <- as.data.frame(fread(paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/mQTLs/Outputs/", i)))
  split <- strsplit(tmp$SNP, "_")
  tmp$MarkerName <- paste(sapply(split, "[[", 1), sapply(split, "[[", 2), sep = ":")
  stats[[i]] <- tmp 
  print(i)
} 


anno = as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Cohort_Summary_Stats/EasyQC/Allele_Freq_and_Mapping_Info/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.rsid_map.gz"))
anno$MarkerName = paste0(anno$chr, ":", anno$pos)


for(i in files){
  stats[[i]]$rsid = anno[match(stats[[i]]$MarkerName, anno$MarkerName), "rsid"]
  write.table(data.frame(stats[[i]]), 
              file = paste0("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/mQTLs/Outputs/", i, ".rsid"), sep='\t', quote=F, row.names=F)
  print(i)
}
