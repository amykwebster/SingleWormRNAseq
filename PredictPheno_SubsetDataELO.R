#start with cpms_phenotype_merge df
library(tidyverse)
load("cpms_phenotype_merge.Rda")
head(cpms_phenotype_merge)
#try to write a model that splits the data in half, estimates top 10 genes from half the data, then estimate ELO on the other half. What is the R2 of that correlation?
brood_genes2<-cpms_phenotype_merge%>%
  dplyr::select(genes,value,ELO_hrs,variable)%>%
  pivot_wider(names_from = "genes",values_from = "value")%>%
  dplyr::select(-variable)


v <- as.vector(c(rep(TRUE,90),rep(FALSE,90))) #create 90 TRUE, 90 FALSE
performance_subset<-c()
for (k in 1:500){ # how many times to split the data in half different ways
  set.seed(k) #make it so random subsetting is reproducible
  ind <- sample(v) #Sample them randomly. 
  brood_subset1 <- brood_genes2[ind, ] 
  brood_subset2 <- brood_genes2[!ind, ]

  highestR2<-c()
  for (i in 2:7939){
  brood_genes3<- brood_subset1 %>%
    dplyr::select(ELO_hrs, i)
  summary_fullModel3<-summary(lm(ELO_hrs~.,brood_genes3))
  highestR2<-rbind(highestR2,c(colnames(brood_genes3)[2],summary_fullModel3$r.squared))
}

highestR2<-as.data.frame(highestR2)
colnames(highestR2)<-c("WBID","R2")
highestR2$R2<-as.numeric(highestR2$R2)
tail(highestR2[order(highestR2$R2),])
tail(highestR2[order(highestR2$R2),])[6,1]
R2_firstgene<-tail(highestR2[order(highestR2$R2),])[6,2]

Gene1<-c(tail(highestR2[order(highestR2$R2),])[6,1])
  for (j in 1:10){ # how many genes included in the predictive model
    highestR2<-c()
    for (i in 2:7939){ # how many genes to loop through to find the best next gene, 2:7939 for all
      brood_genes3<- brood_subset1 %>%
        dplyr::select(ELO_hrs, Gene1[1:j],i)
      summary_fullModel3<-summary(lm(ELO_hrs~.,brood_genes3))
      highestR2<-rbind(highestR2,c(colnames(brood_genes3)[2+j],summary_fullModel3$r.squared))
    }
    highestR2<-as.data.frame(highestR2)
    highestR2$V2<-as.numeric(highestR2$V2)
    colnames(highestR2)<-c("Gene","FullR2")
    highestR2 <- highestR2[order(-highestR2$FullR2),] 
    Gene1<-c(Gene1,highestR2[1,1])
  }
  brood_genes4<- brood_subset1 %>%
    dplyr::select(ELO_hrs, Gene1[1:10]) #CHANGE TO MATCH J
  summary_fullModel4<-summary(lm(ELO_hrs~.,brood_genes4))
  brood_genes5<- brood_subset2 %>%
    dplyr::select(ELO_hrs, Gene1[1:10]) #CHANGE TO MATCH J (and also change in performance_subset)
  summary_fullModel5<-summary(lm(ELO_hrs~.,brood_genes5))
  performance_subset<-rbind(performance_subset,c(Gene1[1:10],summary_fullModel4$r.squared,summary_fullModel5$r.squared))
}



performance_subset<-as.data.frame(performance_subset)
performance_subset$V11<-as.numeric(performance_subset$V11)
performance_subset$V12<-as.numeric(performance_subset$V12)

write.table(performance_subset,"SplitData_Rsquared_ELO_3.txt",quote = F,sep = "\t")





