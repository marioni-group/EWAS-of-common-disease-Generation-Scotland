
# Ideogram 
## Incident findings 
a=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/incident_full_replication_annotations.csv")
a$new_gene=ifelse(nchar(a$new_gene) > 8, paste0(substr(a$new_gene, 1, 7),"."), a$new_gene)
a$new_gene=ifelse(a$Present_in_Basic_Model==1,paste0(a$new_gene,"+"), a$new_gene)
a$new_gene=ifelse(a$Replicated==1,paste0(a$new_gene,"*"), a$new_gene)

a1=a[,c("Probe","Chr", "bp", "trait", "new_gene")]
names(a1)=c("marker", "chrom", "pos", "phenotype", "annotation")
a1[a1$phenotype %in% "diabetes","phenotype"]="Diabetes (Type 2)"
a1[a1$phenotype %in% "heart_disease","phenotype"]="Ischemic Heart Disease"
a1[a1$phenotype %in% "liver_cirrhosis","phenotype"]="Liver Cirrhosis"
a1[a1$phenotype %in% "ovarian_cancer","phenotype"]="Ovarian Cancer"
a1[a1$phenotype %in% "COPD","phenotype"]="COPD"
a1[a1$phenotype %in% "covid_hospitalisation","phenotype"]="COVID Severity"
a1[a1$annotation %in% "LINC006.+*","pos"] = 35320596-25e6
write.table(a1, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/phenogram_test.txt", sep = "\t", row.names=F,quote=F)

## Prevalent findings 

a=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Cleaned_Results/prevalent_full_replication_annotations.csv")
a$new_gene=ifelse(nchar(a$new_gene) > 8, paste0(substr(a$new_gene, 1, 7),"."), a$new_gene)
a$new_gene=ifelse(a$Present_in_Basic_Model==1,paste0(a$new_gene,"+"), a$new_gene)
a$new_gene=ifelse(a$Replicated==1,paste0(a$new_gene,"*"), a$new_gene)

a1=a[,c("Probe","Chr", "bp", "trait", "new_gene")]
names(a1)=c("marker", "chrom", "pos", "phenotype", "annotation")
a1[a1$phenotype %in% "bowel_cancer","phenotype"]="Colorectal Cancer"
a1[a1$phenotype %in% "breast_cancer","phenotype"]="Breast Cancer"
a1[a1$phenotype %in% "CKD","phenotype"]="Chronic Kidney Disease"
a1[a1$phenotype %in% "diabetes","phenotype"]="Diabetes (Type 2)"
a1[a1$phenotype %in% "heart_disease","phenotype"]="Ischemic Heart Disease"
a1[a1$phenotype %in% "parkinsons","phenotype"]="Parkinson's Disease"
a1[a1$phenotype %in% "prostate_cancer","phenotype"]="Prostate Cancer"
a1[a1$phenotype %in% "lung_cancer","phenotype"]="Lung Cancer"
a1[a1$chr %in% 6 & a1$annotation %in% "BAK1","pos"] <- 33383819-15e5
a1[a1$chr %in% 6 & a1$annotation %in% "CUL7","pos"] <- 43004941+15e5
a1[a1$annotation %in% "ERICH1-.","annotation"]="ERICH1-AS1"
write.table(a1, "/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/phenogramprev_test.txt", sep = "\t", row.names=F,quote=F)

# Load requisite libraries 
library(ggplot2)
library(tidyr)
library(cowplot)

element_textbox <- function(...) {
  el <- element_text(...)
  class(el) <- c("element_textbox", class(el))
  el
}

element_grob.element_textbox <- function(element, ...) {
  text_grob <- NextMethod()
  rect_grob <- element_grob(calc_element("strip.background", theme_bw()))
  
  ggplot2:::absoluteGrob(
    grid::gList(
      element_grob(calc_element("strip.background", theme_bw())),
      text_grob
    ),
    height = grid::grobHeight(text_grob), 
    width = grid::unit(1, "npc")
  )
}



## Prevalence 
a=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Results/Results_Final/prevalent_basic.csv")
b=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Results/Results_Final/prevalent_wbcs.csv")
c=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Results/Results_Final/prevalent_full.csv")

a1=as.data.frame(table(a$trait))
names(a1)=c("trait", "basic")
b1=as.data.frame(table(b$trait))
names(b1)=c("trait", "wbcs")
b1=b1[which(!b1$trait %in% "osteo_arthritis"),]
c1=as.data.frame(table(c$trait))
names(c1)=c("trait", "full")
phens=list.files("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Prevalent_Phenotypes/", ".")
phens=phens[which(!phens %in% c("asthma.phen", "parkinsons_family.phen"))]
phens=gsub(".phen", "", phens)
phens1=as.data.frame(phens)
names(phens1)[1]="trait"
phens1$basic=0
phens1$wbcs=0
phens1$full=0
phens1[which(phens1$trait %in% a1$trait),"basic"]<-a1$basic
phens1[which(phens1$trait %in% b1$trait),"wbcs"]<-b1$wbcs
phens1[which(phens1$trait %in% c1$trait),"full"]<-c1$full

test=as.data.frame(gather(phens1, condition, measurement, basic:wbcs:full, factor_key=TRUE))
test=test[!test$trait %in% c("CKD", "pain", "COPD", "diabetes", "heart_disease"),]
test=test[order(test$trait),]
test$paired=rep(1:length(unique(test$trait)),each=3)

x1=ggplot(data=test,aes(x= measurement, y= reorder(factor(trait),measurement))) +
  geom_line(aes(group = paired),color="grey")+
  geom_point(aes(color=factor(condition)), size=4.5) +
  labs(y="Trait")+labs(x="No. of Associations")+ scale_x_continuous(limits=c(0,18),breaks=c(0,2,4,6,8,10,12,14,16,18,20)) +
  theme_classic(24)+
  theme(legend.position="none") +
  scale_color_brewer(palette="Accent", direction=-1) + theme(text = element_text(size=13),axis.text.x = element_text(angle = 45)) + ggtitle("No. of Associations <50") +theme(plot.title = element_text(hjust = 0.5, size = 11))

test=as.data.frame(gather(phens1, condition, measurement, basic:wbcs:full, factor_key=TRUE))
test=test[test$trait %in% c("pain", "COPD", "heart_disease"),]
test=test[order(test$trait),]
test$paired=rep(1:length(unique(test$trait)),each=3)

x2=ggplot(data=test,aes(x= measurement, y= reorder(factor(trait),measurement))) +
  geom_line(aes(group = paired),color="grey")+
  geom_point(aes(color=factor(condition)), size=4.5) +
  labs(y="")+labs(x="No. of Associations")+ scale_x_continuous(limits=c(0,400),breaks=c(0,50,100,150,200,250,300,350,400)) +
  theme_classic(24)+
  theme(legend.position="none") +
  scale_color_brewer(palette="Accent", direction=-1) + theme(text = element_text(size=13),axis.text.x = element_text(angle = 45)) + ggtitle("No. of Associations between 50 and 500") +theme(plot.title = element_text(hjust = 0.5, size = 12.5))


test=as.data.frame(gather(phens1, condition, measurement, basic:wbcs:full, factor_key=TRUE))
test=test[test$trait %in% c("diabetes", "CKD"),]
test=test[order(test$trait),]
test$paired=rep(1:length(unique(test$trait)),each=3)

x3=ggplot(data=test,aes(x= measurement, y= reorder(factor(trait),measurement))) +
  geom_line(aes(group = paired),color="grey")+
  geom_point(aes(color=factor(condition)), size=4.5) +
  labs(y="")+labs(x="No. of Associations")+ scale_x_continuous(limits=c(0,2500),breaks=c(0,250,500,750,1000,1250,1500,1750,2000,2250,2500)) +
  theme_classic(24)+
  theme(legend.position="none") +
scale_color_brewer(palette="Accent", direction=-1) + theme(text = element_text(size=13),axis.text.x = element_text(angle = 45)) +theme(plot.title = element_text(hjust = 0.5, size = 12.5)) 

x4=ggplot(data=test,aes(x= measurement, y= reorder(factor(trait),measurement))) +
  geom_line(aes(group = paired),color="grey")+
  geom_point(aes(color=factor(condition)), size=4.5) +
  labs(y="Trait")+labs(y="Trait")+ scale_x_continuous(limits=c(0,2500),breaks=c(0,250,500,750,1000,1250,1500,1750,2000,2250,2500)) +
  theme_classic(24)+labs("Model",labels=c("Basic","WBC-adjusted", "Fully-adjusted"))+
  theme(legend.position="right") +
  scale_color_brewer(palette="Accent", direction=-1) + theme(text = element_text(size=13),axis.text.x = element_text(angle = 45)) + ggtitle("No. of Associations between >500") +theme(plot.title = element_text(hjust = 0.5))


legend <- get_legend(x4)

## Assemble all components of plot together

plot = ggdraw(plot_grid(x1, plot_grid(x2, x3, nrow = 2, rel_widths = c(1/4,1/4)), legend, ncol = 3, rel_widths = c(0.75,1/2,0.2)))
png("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/prevalence_dumbbell.tiff", units="in", width=14, height=8, res=300)
plot
dev.off()



## Incident 
a=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Results/Results_Final/incident_basic.csv")
b=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Results/Results_Final/incident_wbcs.csv")
c=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Results/Results_Final/incident_full.csv")

a1=as.data.frame(table(a$trait))
names(a1)=c("trait", "basic")
b1=as.data.frame(table(b$trait))
names(b1)=c("trait", "wbcs")
b1=b1[which(!b1$trait %in% "parkinsons"),]
c1=as.data.frame(table(c$trait))
names(c1)=c("trait", "full")
phens=list.files("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/Incident_Phenotypes/", ".")
phens=phens[which(!phens %in% c("10year","5year", "parkinsons.phen"))]
phens=gsub(".phen", "", phens)
phens1=as.data.frame(phens)
names(phens1)[1]="trait"
phens1$basic=0
phens1$wbcs=0
phens1$full=0
phens1[which(phens1$trait %in% a1$trait),"basic"]<-a1$basic
phens1[which(phens1$trait %in% b1$trait),"wbcs"]<-b1$wbcs
phens1[which(phens1$trait %in% c1$trait),"full"]<-c1$full

test=as.data.frame(gather(phens1, condition, measurement, basic:wbcs:full, factor_key=TRUE))
test=test[!test$trait %in% c("COPD", "diabetes", "lung_cancer"),]
test=test[order(test$trait),]
test$paired=rep(1:length(unique(test$trait)),each=3)

x1=ggplot(data=test,aes(x= measurement, y= reorder(factor(trait),measurement))) +
  geom_line(aes(group = paired),color="grey")+
  geom_point(aes(color=factor(condition)), size=4.5) +
  labs(y="Trait")+labs(x="No. of Associations")+ scale_x_continuous(limits=c(0,18),breaks=c(0,2,4,6,8,10,12,14,16,18,20)) +
  theme_classic(24)+
  theme(legend.position="none") +
  scale_color_brewer(palette="Accent", direction=-1) + theme(text = element_text(size=13),axis.text.x = element_text(angle = 45)) + ggtitle("No. of Associations <50") +theme(plot.title = element_text(hjust = 0.5, size = 11))

test=as.data.frame(gather(phens1, condition, measurement, basic:wbcs:full, factor_key=TRUE))
test=test[test$trait %in% c("lung_cancer"),]
test=test[order(test$trait),]
test$paired=rep(1:length(unique(test$trait)),each=3)

x2=ggplot(data=test,aes(x= measurement, y= reorder(factor(trait),measurement))) +
  geom_line(aes(group = paired),color="grey")+
  geom_point(aes(color=factor(condition)), size=4.5) +
  labs(y="")+labs(x="No. of Associations") + scale_x_continuous(limits=c(0,60),breaks=c(0,10,20,30,40,50,60)) +

  theme_classic(24)+
  theme(legend.position="none") +
  scale_color_brewer(palette="Accent", direction=-1) + theme(text = element_text(size=13),axis.text.x = element_text(angle = 45)) + ggtitle("No. of Associations between 50 and 500") +theme(plot.title = element_text(hjust = 0.5, size = 12.5))


test=as.data.frame(gather(phens1, condition, measurement, basic:wbcs:full, factor_key=TRUE))
test=test[test$trait %in% c("diabetes", "COPD"),]
test=test[order(test$trait),]
test$paired=rep(1:length(unique(test$trait)),each=3)

x3=ggplot(data=test,aes(x= measurement, y= reorder(factor(trait),measurement))) +
  geom_line(aes(group = paired),color="grey")+
  geom_point(aes(color=factor(condition)), size=4.5) +
  labs(y="")+labs(x="No. of Associations")+ scale_x_continuous(limits=c(0,14000),breaks=c(0,2000,4000,6000,8000,10000,12000,14000)) +
  theme_classic(24)+
  theme(legend.position="none") +
  scale_color_brewer(palette="Accent", direction=-1) + theme(text = element_text(size=13),axis.text.x = element_text(angle = 45)) +theme(plot.title = element_text(hjust = 0.5, size = 12.5)) 

x4=ggplot(data=test,aes(x= measurement, y= reorder(factor(trait),measurement))) +
  geom_line(aes(group = paired),color="grey")+
  geom_point(aes(color=factor(condition)), size=4.5) +
  labs(y="Trait")+labs(y="Trait")+ scale_x_continuous(limits=c(0,2500),breaks=c(0,250,500,750,1000,1250,1500,1750,2000,2250,2500)) +
  theme_classic(24)+labs("Model",labels=c("Basic","WBC-adjusted", "Fully-adjusted"))+
  theme(legend.position="right") +
  scale_color_brewer(palette="Accent", direction=-1) + theme(text = element_text(size=13),axis.text.x = element_text(angle = 45)) + ggtitle("No. of Associations between >500") +theme(plot.title = element_text(hjust = 0.5))


legend <- get_legend(x4)

## Assemble all components of plot together

plot = ggdraw(plot_grid(x1, plot_grid(x2, x3, nrow = 2, rel_widths = c(1/4,1/4)), legend, ncol = 3, rel_widths = c(0.75,1/2,0.2)))
png("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/incident_dumbbell.tiff", units="in", width=14, height=8, res=300)
plot
dev.off()