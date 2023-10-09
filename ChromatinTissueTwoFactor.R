library(tidyverse)
library(gridExtra)
library(gplots)
setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/combined_analysis/")
Ahringer_ref_table<-read.csv("Ahringer.csv",header = T)
head(Ahringer_ref_table)
length(unique(Ahringer_ref_table$Gene)) #13340

ActGenes<-Ahringer_ref_table$Gene[which(Ahringer_ref_table$Chromatin_domain=="A")]
RegGenes<-Ahringer_ref_table$Gene[which(Ahringer_ref_table$Chromatin_domain=="R")]

GermlineGenes<-Ahringer_ref_table$Gene[grep("Germline", Ahringer_ref_table$Tissue_specificity)]
SomaticClasses<-c("Neurons","Muscle","Intest","Hypod","Soma")
SomaticGenes<-c()
for(i in 1:5){

SomaticGenes<-c(SomaticGenes,Ahringer_ref_table$Gene[grep(SomaticClasses[i],Ahringer_ref_table$Tissue_specificity)])
}


UbiqGenes<-Ahringer_ref_table$Gene[grep("Ubiq",Ahringer_ref_table$Tissue_specificity)]

brood<-read.table("BroodGenes_avgCPMS.csv", sep=",", header=T, stringsAsFactors=F)
head(brood)
venn(list(brood$symbol,brood$sequence,unique(Ahringer_ref_table$Gene)))
4560+2951 #7511


ActB<-rep(0, length=nrow(brood))
RegB<-rep(0, length=nrow(brood))
GermlineB<-rep(0, length=nrow(brood))
SomaticB<-rep(0, length=nrow(brood))
UbiqB<-rep(0, length=nrow(brood))


ActB[which(brood$symbol%in%ActGenes==T|brood$sequence%in%ActGenes==T)]<-1
RegB[which(brood$symbol%in%RegGenes==T|brood$sequence%in%RegGenes==T)]<-2

ActBRegB<-data.frame(ActB,RegB)
head(ActBRegB)
ActBRegB$sumActBRegB<-rowSums(ActBRegB)

ActBRegB%>%
  filter(sumActBRegB>2)%>%
  summarize(totalgenesfiltered=n()) #63

ActBRegB%>%
  filter(sumActBRegB==1)%>%
  summarize(totalgenesfiltered=n()) #4420

ActBRegB%>%
  filter(sumActBRegB==2)%>%
  summarize(totalgenesfiltered=n()) #1509


GermlineB[which(brood$symbol%in%GermlineGenes==T|brood$sequence%in%GermlineGenes==T)]<-1
SomaticB[which(brood$symbol%in%SomaticGenes==T|brood$sequence%in%SomaticGenes==T)]<-2
#UbiqB[which(brood$symbol%in%UbiqGenes==T|brood$sequence%in%UbiqGenes==T)]<-1

SomaGermline<-data.frame(SomaticB,GermlineB)
head(SomaGermline)
SomaGermline$sumSomaGermline<-rowSums(SomaGermline)

SomaGermline%>%
  filter(sumSomaGermline>2)%>%
  summarize(totalgenesfiltered=n()) #528

SomaGermline%>%
  filter(sumSomaGermline==1)%>%
  summarize(totalgenesfiltered=n()) #711

SomaGermline%>%
  filter(sumSomaGermline==2)%>%
  summarize(totalgenesfiltered=n()) #3223

TotalChromatinTissue1<-cbind(brood,ActBRegB)
TotalChromatinTissue<-cbind(TotalChromatinTissue1,SomaGermline)

TotalChromatinTissueSubset<-subset(TotalChromatinTissue,(sumActBRegB==1|sumActBRegB==2)&(sumSomaGermline==1|sumSomaGermline==2))


TotalChromatinTissueSubset$sumSomaGermline <- ifelse(TotalChromatinTissueSubset$sumSomaGermline == 1, "germline", "soma")
TotalChromatinTissueSubset$sumActBRegB <- ifelse(TotalChromatinTissueSubset$sumActBRegB == 1, "active", "regulated")
TotalChromatinTissueSubset_meanError<-TotalChromatinTissueSubset%>%
  group_by(sumSomaGermline,sumActBRegB)%>%
  summarize(meaneffect=mean(CoeffCPMNorm),sdeffect=sd(CoeffCPMNorm))


InteractionEffect<-ggplot(TotalChromatinTissueSubset,aes(x=sumSomaGermline,y=CoeffCPMNorm,color=sumActBRegB))+
  geom_point(alpha=0.1,position=position_jitterdodge(dodge.width = 0.9,jitter.width = 0.2))+
  geom_violin(alpha=0.5)+theme_classic(base_size = 12)+theme(aspect.ratio = 1)+
  geom_point(data=TotalChromatinTissueSubset_meanError,aes(x=sumSomaGermline,y=meaneffect,color=sumActBRegB),position = position_dodge(width = 0.9))+
  labs(x="Tissue",y="Gene effect size on brood")+scale_color_manual(values=c("gray","red"))

InteractionEffect

TotalTissue<-subset(TotalChromatinTissue,sumSomaGermline==1|sumSomaGermline==2)
TotalTissue$sumSomaGermline <- ifelse(TotalTissue$sumSomaGermline == 1, "germline", "soma")
TissueEffect<-ggplot(TotalTissue,aes(x=sumSomaGermline,y=CoeffCPMNorm))+
  geom_jitter(alpha=0.1,width=0.2)+
  geom_violin(alpha=0.5)+theme_classic(base_size = 12)+theme(aspect.ratio = 1.3)+
  labs(x="Tissue",y="Gene effect size on brood")

TotalChromDomain<-subset(TotalChromatinTissue,sumActBRegB==1|sumActBRegB==2)
TotalChromDomain$sumActBRegB <- ifelse(TotalChromDomain$sumActBRegB == 1, "active", "regulated")
ChromDomainEffect<-ggplot(TotalChromDomain,aes(x=sumActBRegB,y=CoeffCPMNorm))+
  geom_jitter(alpha=0.1,width=0.2)+
  geom_violin(alpha=0.5)+theme_classic(base_size = 12)+theme(aspect.ratio = 1.3)+
  labs(x="Chromatin domain",y="Gene effect size on brood")


GraphVenn<-venn(list(UniqueChromDomain=TotalChromDomain$genes,UniqueTissue=TotalTissue$genes))
graphvenn2<-plot(GraphVenn)


summary(glm(CoeffCPMNorm~sumSomaGermline*sumActBRegB,data=TotalChromatinTissueSubset))
summary(aov(CoeffCPMNorm~sumSomaGermline*sumActBRegB,data=TotalChromatinTissueSubset))#same thing

TukeyHSD(aov(CoeffCPMNorm ~ sumSomaGermline*sumActBRegB,data=TotalChromatinTissueSubset), conf.level=0.95)
?TukeyHSD

SomaActive<-subset(TotalChromatinTissueSubset,sumActBRegB=="active" & sumSomaGermline=="soma")
SomaRegulated<-subset(TotalChromatinTissueSubset,sumActBRegB=="regulated" & sumSomaGermline=="soma")
t.test(SomaActive$CoeffCPMNorm, mu = 0, alternative = "two.sided")
t.test(SomaRegulated$CoeffCPMNorm, mu = 0, alternative = "two.sided")


summary(glm(CoeffCPMNorm~sumSomaGermline,data=TotalChromatinTissueSubset))

summary(glm(CoeffCPMNorm~sumActBRegB,data=TotalChromatinTissueSubset))

summary(glm(CoeffCPMNorm~sumSomaGermline,data=TotalTissue))

summary(glm(CoeffCPMNorm~sumActBRegB,data=TotalChromDomain))

TotalChromatinTissueSubset%>%
  group_by(sumSomaGermline,sumActBRegB)%>%
  summarize(total_in_each=n())

TotalChromatinTissue%>%
  group_by(sumActBRegB)%>%
  summarize(total_in_each=n())%>%
  filter(sumActBRegB==3) #63 genes 

TotalChromatinTissue%>%
  group_by(sumSomaGermline)%>%
  summarize(total_in_each=n())%>%
  filter(sumSomaGermline==3) #528

TotalChromatinTissue%>%
  group_by(sumActBRegB)%>%
  summarize(total_in_each=n())%>%
  filter(sumActBRegB==0) #2832 genes 

TotalChromatinTissue%>%
  group_by(sumSomaGermline)%>%
  summarize(total_in_each=n())%>%
  filter(sumSomaGermline==0)

528/(2833+1101+528) #0.118 prop of genes with both soma and germline

grid.arrange(ChromDomainEffect,TissueEffect,InteractionEffect,ncol=2)



