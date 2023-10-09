#start with cpms_phenotype_merge df
load("cpms_phenotype_merge.Rda")
head(cpms_phenotype_merge)
#try to write a model that splits the data in half, estimates top 5 genes from half the data, then estimate brood on the other half. What is the R2 of that correlation?
brood_genes2<-cpms_phenotype_merge%>%
  dplyr::select(genes,value,brood_size,variable)%>%
  pivot_wider(names_from = "genes",values_from = "value")%>%
  dplyr::select(-variable)


v <- as.vector(c(rep(TRUE,90),rep(FALSE,90))) #create 90 TRUE, 90 FALSE
performance_subset<-c()
for (k in 1:50){ # how many times to split the data in half different ways
  set.seed(k) #make it so random subsetting is reproducible
  ind <- sample(v) #Sample them randomly. 
  brood_subset1 <- brood_genes2[ind, ] 
  brood_subset2 <- brood_genes2[!ind, ]
  Gene1<-c("WBGene00000609") #starting gene (col-20)
  for (j in 1:10){ # how many genes included in the predictive model
    highestR2<-c()
    for (i in 2:7939){ # how many genes to loop through to find the best next gene, 2:7939 for all
      brood_genes3<- brood_subset1 %>%
        dplyr::select(brood_size, Gene1[1:j],i)
      summary_fullModel3<-summary(lm(brood_size~.,brood_genes3))
      highestR2<-rbind(highestR2,c(colnames(brood_genes3)[2+j],summary_fullModel3$r.squared))
    }
    highestR2<-as.data.frame(highestR2)
    highestR2$V2<-as.numeric(highestR2$V2)
    colnames(highestR2)<-c("Gene","FullR2")
    highestR2 <- highestR2[order(-highestR2$FullR2),] 
    Gene1<-c(Gene1,highestR2[1,1])
  }
  brood_genes4<- brood_subset1 %>%
    dplyr::select(brood_size, Gene1[1:10]) #CHANGE TO MATCH J
  summary_fullModel4<-summary(lm(brood_size~.,brood_genes4))
  brood_genes5<- brood_subset2 %>%
    dplyr::select(brood_size, Gene1[1:10]) #CHANGE TO MATCH J (and also change in performance_subset)
  summary_fullModel5<-summary(lm(brood_size~.,brood_genes5))
  performance_subset<-rbind(performance_subset,c(Gene1[1:10],summary_fullModel4$r.squared,summary_fullModel5$r.squared))
}


performance_subset
performance_subset<-as.data.frame(performance_subset)
head(performance_subset)
performance_subset$V11<-as.numeric(performance_subset$V11)
performance_subset$V12<-as.numeric(performance_subset$V12)

write.table(performance_subset,"SplitData_Rsquared.txt",quote = F,sep = "\t")


performance_subset<-read.csv("SplitData_Rsquared_Brood_2.csv",header = T)
performance_subset<-read.csv("SplitData_Rsquared_ELO_2.csv",header = T)

head(performance_subset)
mean(performance_subset$V11)
median(performance_subset$V11)
range(performance_subset$V11)

mean(performance_subset$V12)
median(performance_subset$V12)
range(performance_subset$V12)
performance_subset$Iteration<-c(1:50)

PerformanceMelt<-pivot_longer(performance_subset,2:11,names_to = "GeneNumber",values_to = "GeneName")
print(PerformanceMelt,n=31)
PerformanceMelt
colnames(PerformanceMelt)[2:3]<-c("Selected90","Other90")

NumberOfGenes<-PerformanceMelt%>%
  group_by(GeneName)%>%
  summarize(NumberOfGeneRepeats=n())%>%
  arrange(-NumberOfGeneRepeats)%>%
  filter(NumberOfGeneRepeats>2)

setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/singleworm_96samples/")
WS273_geneNames<-read.csv("WS273_geneNames.csv",header = T) #this has the extra information beyond 
setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/combined_analysis/")
NumberOfGenesMerge<-merge(NumberOfGenes,WS273_geneNames,by.x = "GeneName",by.y ="WB_id" )
NumberOfGenesMerge$SymbolSequence<-paste(NumberOfGenesMerge$symbol,NumberOfGenesMerge$sequence,sep = " ")
NumberOfGenesMerge
NumberOfGenesMerge <- NumberOfGenesMerge[order(-NumberOfGenesMerge$NumberOfGeneRepeats),] 
NumberOfGenesMerge$SymbolSequence<-factor(NumberOfGenesMerge$SymbolSequence,levels = NumberOfGenesMerge$SymbolSequence)
TopGenes_Prediction<-ggplot(NumberOfGenesMerge,aes(x=SymbolSequence,y=NumberOfGeneRepeats/200))+
  geom_bar(stat="identity")+theme_bw(base_size = 8)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),aspect.ratio = 1)+
  labs(x="Gene",y="Proportion of iterations with gene selected")+
  ggtitle("Freq. of genes selected")
TopGenes_Prediction
PredictiveValue<-subset(PerformanceMelt,GeneName=="WBGene00000609")
head(PredictiveValueMelt)
PredictiveValueMelt<-pivot_longer(PredictiveValue,1:2,names_to = "SampleSet",values_to = "TotalRSquared")
PredictiveValueMelt$SampleSet<-factor(PredictiveValueMelt$SampleSet,levels = c("Selected90","Other90"))
R2inSplitData<-ggplot(PredictiveValueMelt,aes(SampleSet,TotalRSquared))+
  geom_boxplot()+geom_jitter()+ylim(0,1)+
  theme_bw(base_size = 8)+theme(aspect.ratio = 1)+
  labs(x="Set of worms",y="Total R squared of 10 genes from selected set")+
  ggtitle("R2 of top 10 genes in split data (2 groups of 90)")

gridExtra::grid.arrange(TotalR2_plot,R2inSplitData,TopGenes_Prediction,ncol=3)

head(performance_subset)

