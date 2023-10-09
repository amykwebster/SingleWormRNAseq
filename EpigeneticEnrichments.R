setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/combined_analysis/")
library(tidyverse)
library(gridExtra)
#see if brood genes are more enriched for H3K27me3 than by chance

Ahringer<-read.csv("Ahringer.csv",header = T)

Chromatin_Ahringer<-Ahringer%>%
  filter(Chromatin_domain!="border")%>%
  filter(Gene!=".")%>%
  filter(Chromatin_domain!=".")%>%
  distinct(Gene,Chromatin_domain,.keep_all = TRUE)%>%
  select(Gene,Chromatin_domain)

#make sure each gene is just there once
Chromatin_Ahringer2=Chromatin_Ahringer%>%
  group_by(Gene)%>%
  mutate(number=n())%>%
  filter(number==1)

load("cpms_phenotype_merge.Rda")

GeneSet<-cpms_phenotype_merge%>%
  filter(variable=="102")

setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/singleworm_96samples/")
WS273_geneNames<-read.csv("WS273_geneNames.csv",header = T)
setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/combined_analysis/")

GeneSet_WS273<-merge(GeneSet,WS273_geneNames,by.x = "genes",by.y = "WB_id")

ChromAhringer_sequence<-merge(Chromatin_Ahringer2,GeneSet_WS273,by.x = "Gene",by.y = "sequence")
ChromAhringer_sequence_WBID<-ChromAhringer_sequence[,c(2,4)]
ChromAhringer_symbol<-merge(Chromatin_Ahringer2,GeneSet_WS273,by.x = "Gene",by.y = "symbol")
ChromAhringer_symbol_WBID<-ChromAhringer_symbol[,c(2,4)]
ChromAhringer_all_WBID<-rbind(ChromAhringer_sequence_WBID,ChromAhringer_symbol_WBID)
subset(ChromAhringer_symbol,genes=="WBGene00011558")
subset(ChromAhringer_sequence,genes=="WBGene00011558")

ChromAhringer_all_WBID2<-ChromAhringer_all_WBID%>%
  group_by(genes)%>%
  mutate(sumgene=n())%>%
  filter(sumgene==1)

ChromAhringer_all_WBID2%>%
  group_by(Chromatin_domain)%>%
  summarize(sum_domain=n())%>%
  ungroup()%>%
  mutate(total_genes=sum(sum_domain))%>%
  mutate(prop_domain=sum_domain/total_genes)

#A domain 0.793, R domain 0.207 in background genes that are in Ahringer and my filtered set of almost 8k
#5435 in both

head(ChromAhringer_all_WBID2)
SplitDataBroodAll<-read.csv("SplitData_Rsquared_Brood_2.csv",header = T)
head(SplitDataBroodAll)
SplitDataBroodAllmelt<-SplitDataBroodAll%>%
  select(-V11,-V12)%>%
  pivot_longer(!Iteration,names_to = "GeneOrder",values_to = "Gene")

SplitDataBroodChromAhr<-merge(SplitDataBroodAllmelt,ChromAhringer_all_WBID2,by.x = "Gene",by.y = "genes")
head(SplitDataBroodChromAhr)

SplitDataBroodChromAhr_sum<-SplitDataBroodChromAhr%>%
  group_by(Iteration,Chromatin_domain)%>%
  summarize(sum_domain=n())%>%
  ungroup()%>%
  group_by(Iteration)%>%
  mutate(total_sum=sum(sum_domain))%>%
  mutate(prop_domain=sum_domain/total_sum)

median(SplitDataBroodChromAhr_sum$total_sum) #7
range(SplitDataBroodChromAhr_sum$total_sum) #3 to 10

ggplot(SplitDataBroodChromAhr_sum,aes(x=prop_domain))+
  geom_density(aes(color=Chromatin_domain))+
  geom_vline(xintercept = 0.207)

SplitDataBroodChromAhr_Ronly<-subset(SplitDataBroodChromAhr_sum,Chromatin_domain=="R")
head(SplitDataBroodChromAhr_Ronly)
SplitDataBroodChromAhr_Ronly<-SplitDataBroodChromAhr_Ronly[,c(2,5)]

#ELO
SplitDataELOAll<-read.csv("SplitData_Rsquared_ELO_2.csv",header = T)
head(SplitDataELOAll)
SplitDataELOAllmelt<-SplitDataELOAll%>%
  select(-V11,-V12)%>%
  pivot_longer(!Iteration,names_to = "GeneOrder",values_to = "Gene")

SplitDataELOChromAhr<-merge(SplitDataELOAllmelt,ChromAhringer_all_WBID2,by.x = "Gene",by.y = "genes")
head(SplitDataELOChromAhr)

SplitDataELOChromAhr_sum<-SplitDataELOChromAhr%>%
  group_by(Iteration,Chromatin_domain)%>%
  summarize(sum_domain=n())%>%
  ungroup()%>%
  group_by(Iteration)%>%
  mutate(total_sum=sum(sum_domain))%>%
  mutate(prop_domain=sum_domain/total_sum)

median(SplitDataELOChromAhr_sum$total_sum) #7
range(SplitDataELOChromAhr_sum$total_sum) #3 to 10

ggplot(SplitDataELOChromAhr_sum,aes(x=prop_domain))+
  geom_density(aes(color=Chromatin_domain))+
  geom_vline(xintercept = 0.207)

SplitDataELOChromAhr_Ronly<-subset(SplitDataELOChromAhr_sum,Chromatin_domain=="R")
head(SplitDataELOChromAhr_Ronly)
SplitDataELOChromAhr_Ronly<-SplitDataELOChromAhr_Ronly[,c(2,5)]



#create background set
v <- as.vector(c(rep(TRUE,10),rep(FALSE,7928))) #select groups of 10 from background
GeneSetSampling<-c()
for (k in 1:500){ # how many times to split the data in half different ways
  set.seed(k) #make it so random subsetting is reproducible
  ind <- sample(v) #Sample them randomly. 
  GeneSetSubset <- GeneSet[ind, ] 
  GeneSetSubset$Iteration<-k
  GeneSetSubset<-GeneSetSubset[,c(2,14)]
  GeneSetSampling<-rbind(GeneSetSampling,GeneSetSubset)
}

SplitDataNull<-merge(GeneSetSampling,ChromAhringer_all_WBID2,by= "genes")
head(SplitDataNull)

SplitDataNull_sum<-SplitDataNull%>%
  group_by(Iteration,Chromatin_domain)%>%
  summarize(sum_domain=n())%>%
  ungroup()%>%
  group_by(Iteration)%>%
  mutate(total_sum=sum(sum_domain))%>%
  mutate(prop_domain=sum_domain/total_sum)

head(SplitDataNull_sum)
median(SplitDataNull_sum$total_sum) #7
range(SplitDataNull_sum$total_sum) #2 to 10

SplitDataNull_sum_Ronly<-subset(SplitDataNull_sum,Chromatin_domain=="R")
SplitDataNull_sum_Ronly<-SplitDataNull_sum_Ronly[,c(2,5)]
SplitDataNull_sum_Ronly$GeneSet<-"NullGroup"
SplitDataBroodChromAhr_Ronly$GeneSet<-"PredictiveGroupBrood"
SplitDataELOChromAhr_Ronly$GeneSet<-"PredictiveGroupELO"
SplitDataExpAndNull<-rbind(SplitDataBroodChromAhr_Ronly,SplitDataNull_sum_Ronly)
SplitDataExpAndNull<-rbind(SplitDataExpAndNull,SplitDataELOChromAhr_Ronly)


ggplot(SplitDataExpAndNull,aes(x=prop_domain))+
  geom_density(aes(color=GeneSet),adjust=1)+facet_grid(.~Chromatin_domain)+
  theme_bw()+labs(x="Proportion H3K27me3 domain",y="Density")+
  theme(aspect.ratio = 1)

SplitDataExpAndNull%>%
  group_by(GeneSet)%>%
  summarize(med_prop=median(prop_domain))

PredictBroodHist<-ggplot(subset(SplitDataExpAndNull,GeneSet!="PredictiveGroupELO"),aes(x=prop_domain,fill=GeneSet))+
  geom_histogram(position="identity",alpha=0.5,binwidth = 0.1)+facet_grid(.~Chromatin_domain)+
  theme_classic()+labs(x="Proportion H3K27me3 domain",y="Number of Iterations")+
  theme(aspect.ratio = 1)+
  geom_vline(xintercept = 0.25,alpha=0.5)+
  geom_vline(xintercept = 0.429,alpha=0.5)+
  scale_fill_manual(values = c("black","red"))

PredictBroodHist

PredictELOHist<-ggplot(subset(SplitDataExpAndNull,GeneSet!="PredictiveGroupBrood"),aes(x=prop_domain,fill=GeneSet))+
  geom_histogram(position = "identity",binwidth = 0.1,alpha=0.5)+facet_grid(.~Chromatin_domain)+
  theme_classic()+labs(x="Proportion H3K27me3 domain",y="Number of Iterations")+
  theme(aspect.ratio = 1)+
  geom_vline(xintercept = 0.25,alpha=0.5)+
  geom_vline(xintercept = 0.333,alpha=0.5)+
  scale_fill_manual(values = c("black","blue"))
PredictELOHist

grid.arrange(PredictBroodHist,PredictELOHist)

#statistics for histograms -- look at CDF plot and do KS test
head(SplitDataExpAndNull)


ggplot()+
  stat_ecdf(data = SplitDataExpAndNull,aes(x=prop_domain,color=GeneSet),alpha=0.5)+ #all genes
  labs(y="Cumulative proportion",x="Prop H3K27me3 domain")+
  scale_color_manual(values=c("black","red","blue"))+
  theme_classic(base_size = 10)+theme(aspect.ratio = 1)

EarlyBrood_SplitData<-subset(SplitDataExpAndNull,GeneSet=="PredictiveGroupBrood")
Null_SplitData<-subset(SplitDataExpAndNull,GeneSet=="NullGroup")
ELO_SplitData<-subset(SplitDataExpAndNull,GeneSet=="PredictiveGroupELO")


ks.test(EarlyBrood_SplitData$prop_domain,Null_SplitData$prop_domain,alternative="two.sided") #Kolmogorov Smirnov test

ks.test(ELO_SplitData$prop_domain,Null_SplitData$prop_domain,alternative="two.sided") #Kolmogorov Smirnov test


