

### Prepare GWAS summary data for causal inference analyses #####

# Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/GWAS/")

# Load requisite libraries 
library(data.table)

######################
######## AD ##########
######################

# Read in file 
tmp=as.data.frame(fread("ad.txt"))
# Extract required columns
tmp1=tmp[,c("rsid","chr", "pos", "a1", "a2", "AF", "beta", "se", "P", "Ntotal")]
# Change to standardised format 
names(tmp1)=c("SNP","CHR", "POS","A1", "A2", "FREQ", "BETA", "SE", "P", "N")
# Save out final file 
fwrite(tmp1, "Cleaned/ad.txt",row.names=F)


############################
######## Diabetes ##########
############################

# Read in file 
tmp=as.data.frame(fread("diabetes.txt"))
# Extract required columns 
tmp1=tmp[,c("rsID", "effect_allele", "effect_allele_frequency", "Fixed-effects_beta", "Fixed-effects_SE", "Fixed-effects_p-value")]
# Add in missing columns 
tmp1$N=251739.50
tmp1$A2=NA
tmp1$CHR=NA
tmp1$POS=NA
# Change to standardised format 
names(tmp1)=c("SNP","A1","FREQ","BETA", "SE", "P", "N", "A2", "CHR", "POS")
# Rearrange columns to correct format 
tmp1=tmp1[,c("SNP", "CHR", "POS", "A1", "A2", "FREQ", "BETA", "SE", "P", "N")]
## Additional cleaning 
tmp1$A1=toupper(tmp1$A1)
# Only keep rsids without NA 
tmp1=tmp1[!is.na(tmp1$SNP),]
# Save out file 
fwrite(tmp1, "Cleaned/diabetes.txt",row.names=F)


############################
######## IBD ###############
############################

# Read in file 
tmp=as.data.frame(fread("ibd.txt"))
# Extract required columns 
tmp1=tmp[,c("rsid", "Allele1", "Allele2", "Effect", "StdErr", "P.value")]
# Add in missing columns 
tmp1$N=59957
tmp1$FREQ=NA
tmp1$CHR=NA
tmp1$POS=NA
# Change to standardised format 
names(tmp1)=c("SNP","A1","A2","BETA", "SE", "P", "N", "FREQ", "CHR", "POS")
# Rearrange columns to correct format 
tmp1=tmp1[,c("SNP", "CHR", "POS", "A1", "A2", "FREQ", "BETA", "SE", "P", "N")]
# Additional cleaning 
tmp1$A1=toupper(tmp1$A1)
tmp1$A2=toupper(tmp1$A2)
# Only keep rsids without NA 
tmp1=tmp1[!is.na(tmp1$SNP),]
# Save out file 
fwrite(tmp1, "Cleaned/ibd.txt", row.names=F)



############################
##### Breast Cancer ########
############################

# Read in file 
tmp=as.data.frame(fread("breast_cancer.txt"))
# Extract required columns 
tmp1=tmp[,c("rsid","chr.iCOGs", "Position.iCOGs", "Baseline.Gwas", "Freq.Gwas", "beta.Gwas", "SE.Gwas", "P.value.Gwas")]
# Add in missing columns 
tmp1$N=197755
tmp1$A2=NA
# Change to standardised format 
names(tmp1)=c("SNP","CHR","POS", "A1","FREQ", "BETA", "SE", "P", "N", "A2")
# Rearrange columns to correct format 
tmp1=tmp1[,c("SNP", "CHR", "POS", "A1", "A2", "FREQ", "BETA", "SE", "P", "N")]
# Only keep rsids without NA 
tmp1=tmp1[!is.na(tmp1$SNP),]
# Save out file 
fwrite(tmp1, "Cleaned/breast_cancer.txt", row.names=F)


############################
######### RA ###############
############################

# Read in file 
tmp=as.data.frame(fread("ra_ha.txt"))
# Extract required columns 
tmp1=tmp[,c("rsid", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value", "TotalSampleSize")]
# Add in missing columns 
tmp1$CHR=NA
tmp1$POS=NA
tmp1$A2=NA
# Change to standardised format 
names(tmp1)=c("SNP","A1","A2", "FREQ","BETA", "SE", "P", "N", "CHR", "POS", "A2")
# Rearrange columns to correct format 
tmp1=tmp1[,c("SNP", "CHR", "POS", "A1", "A2", "FREQ", "BETA", "SE", "P", "N")]
# Additional cleaning 
tmp1$A1=toupper(tmp1$A1)
tmp1$A2=toupper(tmp1$A2)
# Only keep rsids without NA 
tmp1=tmp1[!is.na(tmp1$SNP),]
# Save out file 
fwrite(tmp1, "Cleaned/ra.txt", row.names=F)


############################
##### BOWEL CANCER #########
############################

# Read in file 
tmp=as.data.frame(fread("bowel_cancer.txt"))
# Extract required columns 
tmp1=tmp[,c("rsid", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value", "TotalSampleSize")]
# Add in missing columns 
tmp1$CHR=NA
tmp1$POS=NA
tmp1$A2=NA
# Change to standardised format 
names(tmp1)=c("SNP","A1","A2", "FREQ","BETA", "SE", "P", "N", "CHR", "POS", "A2")
# Rearrange columns to correct format 
tmp1=tmp1[,c("SNP", "CHR", "POS", "A1", "A2", "FREQ", "BETA", "SE", "P", "N")]
# Additional cleaning 
tmp1$A1=toupper(tmp1$A1)
tmp1$A2=toupper(tmp1$A2)
# Only keep rsids without NA 
tmp1=tmp1[!is.na(tmp1$SNP),]
# Save out file 
fwrite(tmp1, "Cleaned/bowel_cancer.txt", row.names=F)