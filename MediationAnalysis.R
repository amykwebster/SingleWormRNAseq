library(reshape2)
library(ggplot2)
library(edgeR)
library(nlme)
library(gplots)
library(knitr)
library(DT)
library(tidyverse)
library(RColorBrewer)
library(gplots)
library(psych)
library(matrixStats)


setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/combined_analysis/counts/")
singlewormCountFiles= list.files(pattern= ".txt$",full.names = TRUE)
singlewormCounts=lapply(singlewormCountFiles,read.table, header=T)
setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/singleworm_96samples/")
WS273_geneNames<-read.csv("WS273_geneNames.csv",header = T) #this has the extra information beyond 

setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/combined_analysis/")
singleworm_phenotypes<-read.csv("singleworm_96samples_combined.csv",header = T)
singleworm_phenotypes<-singleworm_phenotypes[,1:11]
head(singleworm_phenotypes)
Sig_brood_genes_448<-read.csv("Significant_brood_combined.csv",header = TRUE)
head(Sig_brood_genes_448)

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
keep_filter<-rowSums(cpm(d2)>1)>=180 #decide how you want to filter the data
d2<-d2[keep_filter,] #filter d2 to only include genes that passed filtering
dim(d2)[1] 
#save(cpm_d2,file="CPM_d2.Rda")
cpm_d2<-cpm(d2,normalized.lib.sizes = TRUE) #make a counts per million object containing normalized CPM
cpm_d2<-as.data.frame(cpm_d2)
cpm_d3<-cpm_d2
cpm_d3$genes<-rownames(cpm_d3)
cpm_d3<-melt(cpm_d3)
cpm_d3$variable<-gsub("_S.*$","",as.character(cpm_d3$variable))
cpms_phenotype_merge<-merge(cpm_d3,singleworm_phenotypes,by.x="variable",by.y="RNA_sample_number")


YA_mediation<-subset(cpms_phenotype_merge,maternal_age=="YA")
OA_mediation<-subset(cpms_phenotype_merge,maternal_age=="OA")
YA_mediation$MA_numeric<-"0"
OA_mediation$MA_numeric<-"1"
mediation_practice1<-rbind(YA_mediation,OA_mediation)
mediation_practice1$MA_numeric<-as.numeric(mediation_practice1$MA_numeric)

lowtemp_mediation<-subset(mediation_practice1,environment=="20C")
hightemp_mediation<-subset(mediation_practice1, environment=="25C_8hr")
lowtemp_mediation$temp_numeric<-"0"
hightemp_mediation$temp_numeric<-"1"
mediation_practice1<-rbind(lowtemp_mediation,hightemp_mediation)
class(mediation_practice1$temp_numeric)
mediation_practice1$temp_numeric<-as.numeric(mediation_practice1$temp_numeric)
dim(mediation_practice1)
1428840/180


library(psych)

Split_CPMs<-split(mediation_practice1,f=mediation_practice1$genes)
lme_brood_mediate<-c()
for (i in 1:length(Split_CPMs)){
  mediation_output<-mediate(brood_size ~ value + (MA_numeric), data = Split_CPMs[[i]], n.iter = 1000,plot = FALSE)
  lme_brood_mediate<-rbind(lme_brood_mediate,c(Split_CPMs[[i]]$genes[1],mediation_output$c[1,1],mediation_output$cprime[2,1],mediation_output$total.reg$prob[2,1],mediation_output$cprime.reg$prob[2,1],mediation_output$boot$ci.ab[1,1],mediation_output$boot$ci.ab[2,1],mediation_output$a.reg$prob[2,1]))
}

lme_brood_mediate<-as.data.frame(lme_brood_mediate)
head(lme_brood_mediate)
colnames(lme_brood_mediate)<-c("Gene","Model_GeneExpOnly","Model_WithMA","Pval_GeneExp","Pval_WithMA","CI_lowbound","CI_highbound","Pval_mediator")
lme_brood_mediate1<-lme_brood_mediate
lme_brood_mediate$Model_GeneExpOnly<-as.numeric(as.character(lme_brood_mediate$Model_GeneExpOnly))
lme_brood_mediate[,3]<-as.numeric(as.character(lme_brood_mediate[,3]))
lme_brood_mediate[,4]<-as.numeric(as.character(lme_brood_mediate[,4]))
lme_brood_mediate[,5]<-as.numeric(as.character(lme_brood_mediate[,5]))
lme_brood_mediate[,6]<-as.numeric(as.character(lme_brood_mediate[,6]))
lme_brood_mediate[,7]<-as.numeric(as.character(lme_brood_mediate[,7]))
lme_brood_mediate[,8]<-as.numeric(as.character(lme_brood_mediate[,8]))

write.table(lme_brood_mediate,"MediationAnalysis_Brood.txt",quote = F,sep = "\t")

lme_brood_mediate<-read.csv("MediationAnalysis_Brood.csv",header = T)

lme_brood_mediate$Prop_GeneExp<- (abs(lme_brood_mediate$Model_WithMA) / abs(lme_brood_mediate$Model_GeneExpOnly))
lme_brood_mediate$Prop_Mediated<- ((abs(lme_brood_mediate$Model_GeneExpOnly) - abs(lme_brood_mediate$Model_WithMA)) / abs(lme_brood_mediate$Model_GeneExpOnly))

head(lme_brood_mediate)
head(BroodFullDE)
SignificantLMEBroodMediate1<-merge(lme_brood_mediate,BroodFullDE,by= "Gene")
head(SignificantLMEBroodMediate1)
SignificantLMEBroodMediate<-subset(SignificantLMEBroodMediate1,PVal_GeneExp< 0.05/8824) #672
head(SignificantLMEBroodMediate)
range(SignificantLMEBroodMediate$Prop_Mediated)

Sig_Only_Mediations<-merge(mediation_practice1,Sig_brood_genes_448,by.x = "genes",by.y = "Gene")
77400/180

SignificantLMEBroodMediate <- SignificantLMEBroodMediate[order(-SignificantLMEBroodMediate$Prop_Mediated),] 
head(SignificantLMEBroodMediate)
SignificantLMEBroodMediate$Gene<-factor(SignificantLMEBroodMediate$Gene,levels = SignificantLMEBroodMediate$Gene)
ggplot(SignificantLMEBroodMediate,aes(x=Gene,y=Prop_GeneExp))+
  geom_point(aes())+coord_flip()+theme_classic()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 1)+
  labs(x="Significant Gene",y="Proportion of gene expression effect independent of maternal age")


ggplot(SignificantLMEBroodMediate,aes(x=Prop_Mediated))+
  geom_histogram()+theme_bw()+labs(x="Prop. maternal effect",y="# of genes")+
  theme(aspect.ratio = 1.3)



subset(SignificantLMEBroodMediate,Prop_GeneExp <0.1)
subset(SignificantLMEBroodMediate,Prop_GeneExp >0.8)

median(SignificantLMEBroodMediate$Prop_GeneExp)
subset(SignificantLMEBroodMediate,Gene=="WBGene00002010")


SigBrood_SigMediation<-subset(SignificantLMEBroodMediate,Pval_mediator>0.1)
tail(SigBrood_SigMediation)
range(SigBrood_SigMediation$Prop_Mediated)
SigBrood_SigMediation <- SigBrood_SigMediation[order(SigBrood_SigMediation$Prop_Mediated),] 

mediation_output<-mediate(brood_size ~ value + (MA_numeric), data = Split_CPMs[[726]], n.iter = 1000)
mediation_output
tail(lme_brood_mediate)
head(lme_brood_mediate)



mediate(brood_size ~ value + (MA_numeric), data = Split_CPMs[[i]], n.iter = 1000)%>% print(short = FALSE)

mediation_output<-mediate(brood_size ~ value + (MA_numeric), data = Split_CPMs[[i]], n.iter = 1000,plot = FALSE)

mediation_output$c[1,1] #total model c
mediation_output$cprime[2,1] #c prime additive model
mediation_output$total.reg$prob[2,1] #p-val for gene exp only
mediation_output$cprime.reg$prob[2,1] #p-val for gene exp when mat age included
mediation_output$boot$ci.ab[1,1] #2.5% confidence of (c-c')
mediation_output$boot$ci.ab[2,1]#97.5% confidence of (c-c')
mediation_output$a.reg$prob[2,1] #p-val of mediator on gene exp


head(Split_CPMs[[1]]$symbol[1])

#Repeat but only for 430 genes (of the 448)
head(Sig_Only_Mediations)
head(Split_CPMs[[1]])
Split_CPMs<-split(Sig_Only_Mediations,f=Sig_Only_Mediations$genes)
lme_brood_mediate1<-c()
for (i in 1:length(Split_CPMs)){
  mediation_output<-mediate(brood_size ~ value + (MA_numeric), data = Split_CPMs[[i]], n.iter = 1000,plot = FALSE)
  lme_brood_mediate1<-rbind(lme_brood_mediate1,c(Split_CPMs[[i]]$genes[1],Split_CPMs[[i]]$symbol[1],Split_CPMs[[i]]$sequence[1],mediation_output$c[1,1],mediation_output$cprime[2,1],mediation_output$total.reg$prob[2,1],mediation_output$cprime.reg$prob[2,1],mediation_output$boot$ci.ab[1,1],mediation_output$boot$ci.ab[2,1],mediation_output$a.reg$prob[2,1],mediation_output$ab[1,1]))
}
mediation_output<-mediate(brood_size ~ value + (MA_numeric), data = Split_CPMs[[430]], n.iter = 1000,plot = FALSE)
mediation_output$c
mediation_output$cprime
mediation_output$ab
mediation_output$all.ab

c(Split_CPMs[[i]]$genes[1],Split_CPMs[[430]]$symbol[1],Split_CPMs[[i]]$sequence[1],mediation_output$c[1,1],mediation_output$cprime[2,1],mediation_output$total.reg$prob[2,1],mediation_output$cprime.reg$prob[2,1],mediation_output$boot$ci.ab[1,1],mediation_output$boot$ci.ab[2,1],mediation_output$a.reg$prob[2,1])

lme_brood_mediate1<-as.data.frame(lme_brood_mediate1)
tail(lme_brood_mediate1)
colnames(lme_brood_mediate1)<-c("Gene","Symbol","Sequence","TotalEffect_GeneExp","DirectEffect_GeneExp","Pval_GeneExp","Pval_WithMA","CI_lowbound","CI_highbound","Pval_mediator","IndirectEffect_GeneExp")
lme_brood_mediate2<-lme_brood_mediate1
lme_brood_mediate1$TotalEffect_GeneExp<-as.numeric(as.character(lme_brood_mediate1$TotalEffect_GeneExp))
lme_brood_mediate1[,5]<-as.numeric(as.character(lme_brood_mediate1[,5]))
lme_brood_mediate1[,6]<-as.numeric(as.character(lme_brood_mediate1[,6]))
lme_brood_mediate1[,7]<-as.numeric(as.character(lme_brood_mediate1[,7]))
lme_brood_mediate1[,8]<-as.numeric(as.character(lme_brood_mediate1[,8]))
lme_brood_mediate1[,9]<-as.numeric(as.character(lme_brood_mediate1[,9]))
lme_brood_mediate1[,10]<-as.numeric(as.character(lme_brood_mediate1[,10]))
lme_brood_mediate1[,11]<-as.numeric(as.character(lme_brood_mediate1[,11]))

0.02898075/0.05835169
0.007040877/0.05835169



write.table(lme_brood_mediate1,"MediationAnalysis_Brood_430of448.txt",quote = F,sep = "\t")

lme_brood_mediate1$Prop_Direct<- (abs(lme_brood_mediate1$DirectEffect_GeneExp) / abs(lme_brood_mediate1$TotalEffect_GeneExp))
lme_brood_mediate1$Prop_Mediated<- (abs(lme_brood_mediate1$IndirectEffect_GeneExp) / abs(lme_brood_mediate1$TotalEffect_GeneExp))
lme_brood_mediate1$Mediate_LowBound<-(abs(lme_brood_mediate1$CI_lowbound) / abs(lme_brood_mediate1$TotalEffect_GeneExp))
lme_brood_mediate1$Mediate_HighBound<-(abs(lme_brood_mediate1$CI_highbound) / abs(lme_brood_mediate1$TotalEffect_GeneExp))
lme_brood_mediate1$Sum<-lme_brood_mediate1$DirectEffect_GeneExp + lme_brood_mediate1$IndirectEffect_GeneExp
lme_brood_mediate1$Difference<-lme_brood_mediate1$Sum - lme_brood_mediate1$TotalEffect_GeneExp

head(lme_brood_mediate1)
range(lme_brood_mediate1$Difference)#yes this works! "Sum" should equal total effect, "Difference" should be close to 0


0.4966566-0.3026888
0.3026888-0.1206628

SignificantLMEBroodMediate1<-subset(lme_brood_mediate1,Pval_GeneExp<0.1) 
head(SignificantLMEBroodMediate1)
range(SignificantLMEBroodMediate1$Prop_Mediated)
subset(SignificantLMEBroodMediate1,Pval_WithMA>0.1)

UPgenes<-subset(SignificantLMEBroodMediate1,TotalEffect_GeneExp>0)
DOWNgenes<-subset(SignificantLMEBroodMediate1,TotalEffect_GeneExp<0)

UPgenes$Direction<-"PositiveAssoc"
DOWNgenes$Direction<-"NegativeAssoc"

SignificantLMEBroodMediate2<-rbind(UPgenes,DOWNgenes)

SignificantLMEBroodMediate2 <- SignificantLMEBroodMediate2[order(-SignificantLMEBroodMediate2$Prop_Mediated),] 
head(SignificantLMEBroodMediate2)
SignificantLMEBroodMediate2$Gene<-factor(SignificantLMEBroodMediate2$Gene,levels = SignificantLMEBroodMediate2$Gene)
MAonBrood<-ggplot(SignificantLMEBroodMediate2,aes(x=Gene,y=Prop_Mediated))+
  facet_grid(.~Direction)+geom_hline(yintercept=0.5)+ylim(-0.2,1.25)+
  #geom_point(aes(x=Gene,y=Mediate_LowBound),color="blue")+
  #geom_point(aes(x=Gene,y=Mediate_HighBound),color="red")+
  geom_errorbar(aes(x=Gene,ymin=Mediate_LowBound,ymax=Mediate_HighBound,color=Direction),width=0.05,size=1,alpha=0.1)+
  geom_point(aes())+
  coord_flip()+theme_classic()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 1)+
  labs(x="Significant Gene",y="Proportion of gene expression effect dependent on maternal age")
MAonBrood

head(SignificantLMEBroodMediate2)
df$x <- paste(df$n,df$s)
SignificantLMEBroodMediate2$SymbolSequence<-paste(SignificantLMEBroodMediate2$Symbol,SignificantLMEBroodMediate2$Sequence)
head(SignificantLMEBroodMediate2)
MostDependent<-subset(SignificantLMEBroodMediate2,Prop_Mediated>0.75)

MostIndependent<-subset(SignificantLMEBroodMediate2,Prop_Mediated<0.25)

MostDependent$SymbolSequence<-factor(MostDependent$SymbolSequence,levels = MostDependent$SymbolSequence)
MostDep<-ggplot(MostDependent,aes(x=SymbolSequence,y=Prop_Mediated))+
  geom_point(aes())+
  #geom_point(aes(x=Gene,y=Mediate_LowBound),color="blue")+
  #geom_point(aes(x=Gene,y=Mediate_HighBound),color="red")+
  geom_errorbar(aes(x=SymbolSequence,ymin=Mediate_LowBound,ymax=Mediate_HighBound,color=Direction),width=0.2,size=1,alpha=0.5)+
  coord_flip()+theme_classic(base_size = 10)+theme(aspect.ratio = 1)+
  labs(x="Significant Gene",y="Proportion of gene expression effect dependent on maternal age")+
  ggtitle("> 75% effect driven by maternal age")

MostIndependent$SymbolSequence<-factor(MostIndependent$SymbolSequence,levels = MostIndependent$SymbolSequence)
MostInd<-ggplot(MostIndependent,aes(x=SymbolSequence,y=Prop_Mediated))+
  geom_point(aes())+
  #geom_point(aes(x=Gene,y=Mediate_LowBound),color="blue")+
  #geom_point(aes(x=Gene,y=Mediate_HighBound),color="red")+
  geom_errorbar(aes(x=SymbolSequence,ymin=Mediate_LowBound,ymax=Mediate_HighBound,color=Direction),width=0.2,size=1,alpha=0.5)+
  coord_flip()+theme_classic(base_size = 10)+theme(aspect.ratio = 1)+
  labs(x="Significant Gene",y="Proportion of gene expression effect dependent on maternal age")+
  ggtitle("< 25% effect driven by maternal age") 

gridExtra::grid.arrange(MostDep,MostInd,ncol=2)

ggplot(subset(SignificantLMEBroodMediate2,Prop_GeneExp<0.25),aes(x=Gene,y=Prop_GeneExp))+
  geom_point(aes(color=Direction))+coord_flip()+theme_classic()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 1)+
  geom_text(aes(label=Symbol),nudge_x = 0.33)+facet_grid(.~Direction)+
  labs(x="Significant Gene",y="Proportion of gene expression effect independent of maternal age")


ggplot(subset(SignificantLMEBroodMediate2,Prop_GeneExp>0.75),aes(x=Gene,y=Prop_GeneExp))+
  geom_point(aes(color=Direction))+coord_flip()+theme_classic()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 1)+
  geom_text(aes(label=Symbol),nudge_x = 0.33,size=4)+facet_grid(.~Direction)+
  labs(x="Significant Gene",y="Proportion of gene expression effect independent of maternal age")


subset(SignificantLMEBroodMediate2,Symbol=="puf-5" | Symbol=="puf-7")
head(cpms_phenotype_merge)
head(MostDependent)
MostDependent_cpmsPhenotype<-merge(cpms_phenotype_merge,MostDependent,by.x = "genes",by.y = "Gene")
head(MostDependent_cpmsPhenotype)
MostIndependent_cpmsPhenotype<-merge(cpms_phenotype_merge,MostIndependent,by.x = "genes",by.y = "Gene")

MostDep_Exp<-ggplot(MostDependent_cpmsPhenotype,aes(x=value,y=brood_size))+
  geom_point(alpha=0.5,size=0.5)+facet_grid(.~SymbolSequence,scales = "free_x")+
  geom_smooth(method = "lm")+
  theme_bw(base_size = 5)+
  ggtitle("Most dependent genes")+
  labs(x="CPMs",y="Early brood size")+
  theme(aspect.ratio = 1)


MostInd_Exp<-ggplot(MostIndependent_cpmsPhenotype,aes(x=value,y=brood_size))+
  geom_point(alpha=0.5,size=0.5)+facet_grid(.~SymbolSequence,scales = "free_x")+
  geom_smooth(method = "lm")+
  theme_bw(base_size = 5)+
  ggtitle("Most independent genes")+
  labs(x="CPMs",y="Early brood size")+
  theme(aspect.ratio = 1)

MostDep_Exp_MA<-ggplot(MostDependent_cpmsPhenotype,aes(x=value,y=brood_size,color=maternal_age))+
  geom_point(alpha=0.5,size=0.5)+facet_grid(.~SymbolSequence,scales = "free_x")+
  geom_smooth(method = "lm")+
  theme_bw(base_size = 5)+
  ggtitle("Most dependent genes")+
  labs(x="CPMs",y="Early brood size")+
  theme(aspect.ratio = 1)

MostInd_Exp_MA<-ggplot(MostIndependent_cpmsPhenotype,aes(x=value,y=brood_size,color=maternal_age))+
  geom_point(alpha=0.5,size=0.5)+facet_grid(.~SymbolSequence,scales = "free_x")+
  geom_smooth(method = "lm")+
  theme_bw(base_size = 5)+
  ggtitle("Most independent genes")+
  labs(x="CPMs",y="Early brood size")+
  theme(aspect.ratio = 1)

head(SignificantLMEBroodMediate2)
gridExtra::grid.arrange(MostDep_Exp,MostDep_Exp_MA,MostInd_Exp,MostInd_Exp_MA,ncol=1)

#Split data into 4 groups
UP_independent<-subset(SignificantLMEBroodMediate2,TotalEffect_GeneExp >0 & Prop_Direct > 0.5) #66
DOWN_independent<-subset(SignificantLMEBroodMediate2,TotalEffect_GeneExp <0 & Prop_Direct > 0.5) #202
UP_dependent<-subset(SignificantLMEBroodMediate2,TotalEffect_GeneExp >0 & Prop_Direct < 0.5) #92
DOWN_dependent<-subset(SignificantLMEBroodMediate2,TotalEffect_GeneExp <0 & Prop_Direct < 0.5) #70

write.table(UP_independent,"UP_independent_66.txt",quote = F,sep = "\t")
write.table(DOWN_independent,"DOWN_independent_202.txt",quote = F,sep = "\t")
write.table(UP_dependent,"UP_dependent_92.txt",quote = F,sep = "\t")
write.table(DOWN_dependent,"DOWN_dependent_70.txt",quote = F,sep = "\t")

66+202+92+70
head(Split_CPMs[[1]])
class(Split_CPMs[[1]]$brood_size)
class(Split_CPMs[[1]]$ELO_hrs)

#make sure they aren't dependent on early-life temperature
Split_CPMs<-split(Sig_Only_Mediations,f=Sig_Only_Mediations$genes)
lme_brood_mediate3<-c()
for (i in 1:length(Split_CPMs)){
  mediation_output<-mediate(ELO_hrs ~ value + (temp_numeric), data = Split_CPMs[[i]], n.iter = 1000,plot = FALSE)
  lme_brood_mediate3<-rbind(lme_brood_mediate3,c(Split_CPMs[[i]]$genes[1],Split_CPMs[[i]]$symbol[1],Split_CPMs[[i]]$sequence[1],mediation_output$c[1,1],mediation_output$cprime[2,1],mediation_output$total.reg$prob[2,1],mediation_output$cprime.reg$prob[2,1],mediation_output$boot$ci.ab[1,1],mediation_output$boot$ci.ab[2,1],mediation_output$a.reg$prob[2,1]))
}
mediation_output
lme_brood_mediate3<-as.data.frame(lme_brood_mediate3)
head(lme_brood_mediate3)
colnames(lme_brood_mediate3)<-c("Gene","Symbol","Sequence","Model_GeneExpOnly","Model_WithTemp","Pval_GeneExp","Pval_WithTemp","CI_lowbound","CI_highbound","Pval_mediator")
lme_brood_mediate4<-lme_brood_mediate3
lme_brood_mediate3$Model_GeneExpOnly<-as.numeric(as.character(lme_brood_mediate3$Model_GeneExpOnly))
lme_brood_mediate3[,5]<-as.numeric(as.character(lme_brood_mediate3[,5]))
lme_brood_mediate3[,6]<-as.numeric(as.character(lme_brood_mediate3[,6]))
lme_brood_mediate3[,7]<-as.numeric(as.character(lme_brood_mediate3[,7]))
lme_brood_mediate3[,8]<-as.numeric(as.character(lme_brood_mediate3[,8]))
lme_brood_mediate3[,9]<-as.numeric(as.character(lme_brood_mediate3[,9]))
lme_brood_mediate3[,10]<-as.numeric(as.character(lme_brood_mediate3[,10]))


lme_brood_mediate3$Prop_GeneExp<- (abs(lme_brood_mediate3$Model_WithTemp) / abs(lme_brood_mediate3$Model_GeneExpOnly))
lme_brood_mediate3$Prop_Mediated<- ((abs(lme_brood_mediate3$Model_GeneExpOnly) - abs(lme_brood_mediate3$Model_WithTemp)) / abs(lme_brood_mediate3$Model_GeneExpOnly))
lme_brood_mediate3$Mediate_LowBound<-((lme_brood_mediate3$CI_lowbound) / (lme_brood_mediate3$Model_GeneExpOnly))
lme_brood_mediate3$Mediate_HighBound<-((lme_brood_mediate3$CI_highbound) / (lme_brood_mediate3$Model_GeneExpOnly))



SignificantLMEBroodMediate3<-subset(lme_brood_mediate3,Pval_GeneExp<0.1) 
head(SignificantLMEBroodMediate3)
range(SignificantLMEBroodMediate3$Prop_Mediated) #no more than 5%

UPgenes3<-subset(SignificantLMEBroodMediate3,Model_GeneExpOnly>0)
DOWNgenes3<-subset(SignificantLMEBroodMediate3,Model_GeneExpOnly<0)

UPgenes3$Direction<-"PositiveAssoc"
DOWNgenes3$Direction<-"NegativeAssoc"

SignificantLMEBroodMediate4<-rbind(UPgenes3,DOWNgenes3)

SignificantLMEBroodMediate4 <- SignificantLMEBroodMediate4[order(-SignificantLMEBroodMediate4$Prop_Mediated),] 
head(SignificantLMEBroodMediate4)
SignificantLMEBroodMediate4$Gene<-factor(SignificantLMEBroodMediate4$Gene,levels = SignificantLMEBroodMediate4$Gene)
TempOnBrood<-ggplot(SignificantLMEBroodMediate4,aes(x=Gene,y=Prop_Mediated))+
  facet_grid(.~Direction)+geom_hline(yintercept=0.5)+ylim(-0.2,1.25)+
  #geom_point(aes(x=Gene,y=Mediate_LowBound),color="blue")+
  #geom_point(aes(x=Gene,y=Mediate_HighBound),color="red")+
  geom_errorbar(aes(x=Gene,ymin=Mediate_LowBound,ymax=Mediate_HighBound,color=Direction),width=0.05,size=1,alpha=0.1)+
  geom_point(aes())+
  coord_flip()+theme_classic()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 1)+
  labs(x="Significant Gene",y="Proportion of gene expression effect dependent on temperature")

gridExtra::grid.arrange(MAonBrood,TempOnBrood,ncol=1)


