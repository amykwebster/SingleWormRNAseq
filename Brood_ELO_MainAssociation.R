library(reshape2)
library(ggplot2)
library(edgeR)
library(nlme)
library(gplots)
library(knitr)
library(DT)
library(tidyverse)
library(RColorBrewer)

setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/ScriptsFilesForPaper/")
singleworm_phenotypes<-read.csv("singleworm_192samples.csv",header = T)
singleworm_phenotypes180worms<-subset(singleworm_phenotypes,RNA_sample_number!="101"& RNA_sample_number!="178"& RNA_sample_number!="179"& RNA_sample_number!="76" & RNA_sample_number!="80" & RNA_sample_number!="85"& RNA_sample_number!="86" & RNA_sample_number!="87" & RNA_sample_number!="94" & RNA_sample_number!="99" & RNA_sample_number!="355" & RNA_sample_number!="359")

setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/combined_analysis/counts/")
singlewormCountFiles= list.files(pattern= ".txt$",full.names = TRUE)
singlewormCounts=lapply(singlewormCountFiles,read.table, header=T)

library1<-singlewormCounts[[1]][,c(1,7)]
library2<-singlewormCounts[[2]][,c(1,7)]
df_counts<-merge(library1,library2,by="Geneid")

for (i in 3:193){
  df_counts<-merge(df_counts,singlewormCounts[[i]][,c(1,7)])
}

df_counts1<-df_counts
names(df_counts)<-gsub("X20220222_singleworm.bams.","",names(df_counts))
names(df_counts)<-gsub("_L001_R1_001.bam","",names(df_counts))
names(df_counts)<-gsub("X20221220_singleworm_batch2.bams.","",names(df_counts))
names(df_counts)<-gsub("_L002_R1_001.bam","",names(df_counts))
rownames(df_counts)<-df_counts$Geneid
df_counts<-df_counts[,c(-1,-194)] #Gets it down to the counts for the 96 samples

setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/singleworm_96samples/")
WS273_geneNames<-read.csv("WS273_geneNames.csv",header = T) #this has the extra information beyond the wormbase IDs (WB_id)

#Cut down to protein-coding genes
counts_merge<-merge(df_counts,WS273_geneNames,by.x = 0,by.y = "WB_id") 
rownames(counts_merge)<-counts_merge$Row.names
df_counts_proteinCoding<-subset(counts_merge,type=="protein_coding_gene") #20095 protein coding genes, lines up with what wormbase says
df_counts2<-df_counts_proteinCoding[,2:193]


df_counts3<-df_counts2
df_counts3$genes<-rownames(df_counts3)
df_counts3<-melt(df_counts3)
df_counts3$variable<-gsub("_S.*$","",as.character(df_counts3$variable))
counts_phenotype_merge<-merge(df_counts3,singleworm_phenotypes,by.x="variable",by.y="RNA_sample_number")
counts_phenotype_merge_1gene<-subset(counts_phenotype_merge,genes=="WBGene00000001")
#Remove 12 libraries that had low sequencing coverage and/or other issues
counts_phenotype_merge_1gene<-subset(counts_phenotype_merge_1gene,variable!="101"& variable!="178"& variable!="179"& variable!="76" & variable!="80" & variable!="85"& variable!="86" & variable!="87" & variable!="94" & variable!="99" & variable!="355" & variable!="359")

#Again remove 10 libraries that had low sequencing coverage from the counts table
df_counts2<-select(df_counts2,-`101_S30`,-`178_S59`,-`179_S60`,-`76_S9`,-`80_S13`,-`85_S17`,-`86_S18`,-`87_S19`,-`94_S24`,-`99_S29`,-`355_S48`,-`359_S96`)

#Using edgeR for differential expression analysis
counts_groups_AGE<-as.character(counts_phenotype_merge_1gene$maternal_age)
counts_groups_TEMP<-as.character(counts_phenotype_merge_1gene$environment)

#Create DGE object in edgeR based on the maternal age environmental condition
d2<-DGEList(counts=df_counts2,group=factor(counts_groups_AGE))
#check the dimensions of the d2 without any filtering
keep_filter<-rowSums(cpm(d2)>10)>=1 #decide how you want to filter the data
#keep_filter<-rowSums(cpm(d2)>1)>=180 #ALTERNATE FILTERING!!!

d2<-d2[keep_filter,] #filter d2 to only include genes that passed filtering
dim(d2)[1] #8824 genes, or 7938 for alternate filtering
cpm_d2<-cpm(d2,normalized.lib.sizes = TRUE) #make a counts per million object containing normalized CPM
cpm_d2<-as.data.frame(cpm_d2)
cpm_d3<-cpm_d2
cpm_d3$genes<-rownames(cpm_d3)
cpm_d3<-melt(cpm_d3)
cpm_d3$variable<-gsub("_S.*$","",as.character(cpm_d3$variable))
cpms_phenotype_merge<-merge(cpm_d3,singleworm_phenotypes,by.x="variable",by.y="RNA_sample_number")

#ANALYSIS: Gene expression variation underlying brood size and egg-laying onset
Split_CPMs<-split(cpms_phenotype_merge,f=cpms_phenotype_merge$genes)
lme_brood<-c()
for (i in 1:length(Split_CPMs)){
  summary_brood_model<-summary(lme(brood_size~ value,random=~1|replicate,data = Split_CPMs[[i]]))
  lme_brood<-rbind(lme_brood,c(Split_CPMs[[i]]$genes[1],summary_brood_model$tTable[2,5],summary_brood_model$tTable[2,1]))
}
lme_brood<-as.data.frame(lme_brood)
colnames(lme_brood)<-c("Gene","PVal_Exp","UpOrDown")
lme_brood$PVal_GeneExp<-as.numeric(as.character(lme_brood$PVal_Exp))
lme_brood$UpDown<-as.numeric(as.character(lme_brood$UpOrDown))
#write.table(lme_brood,"BroodGenes_8824.txt",quote = F,sep = "\t")
lmeBrood_SigExp<-subset(lme_brood,PVal_GeneExp<5.666364e-06)
lmeBrood_SigExp_geneNames<-merge(lmeBrood_SigExp,WS273_geneNames,by.x = "Gene",by.y="WB_id")
UpBrood<-subset(lmeBrood_SigExp_geneNames,UpDown>0)
DownBrood<-subset(lmeBrood_SigExp_geneNames,UpDown<0)
sigBroodTable<-lmeBrood_SigExp_geneNames[,c(-2,-3,-8,-9)]
colnames(sigBroodTable)<-c("Gene","P-value","Magnitude/direction","Symbol","Sequence")

lme_ELO<-c()
for (i in 1:length(Split_CPMs)){
  summary_brood_model<-summary(lme(ELO_hrs ~ value,random=~1|replicate,data = Split_CPMs[[i]]))
  lme_ELO<-rbind(lme_ELO,c(Split_CPMs[[i]]$genes[1],summary_brood_model$tTable[2,5],summary_brood_model$tTable[2,1]))
}

#11 genes
lme_ELO<-as.data.frame(lme_ELO)
colnames(lme_ELO)<-c("Gene","PVal_Exp","UpOrDown")
lme_ELO$PVal_GeneExp<-as.numeric(as.character(lme_ELO$PVal_Exp))
lme_ELO$UpDown<-as.numeric(as.character(lme_ELO$UpOrDown))
lmeELO_SigExp<-subset(lme_ELO,PVal_GeneExp< 5.666364e-06)
lmeELO_SigExp_geneNames<-merge(lmeELO_SigExp,WS273_geneNames,by.x = "Gene",by.y="WB_id")
UpELO<-subset(lmeELO_SigExp_geneNames,UpDown>0)
DownELO<-subset(lmeELO_SigExp_geneNames,UpDown<0)
sigELOTable<-lmeELO_SigExp_geneNames[,c(-2,-3,-8,-9)]
colnames(sigELOTable)<-c("Gene","P-value","Magnitude/direction","Symbol","Sequence")



