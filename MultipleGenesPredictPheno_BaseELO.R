#Next, INDIVIDUAL GENES rather than Principal Components
#Baseline for ELO- non shuffled data
library(tidyverse)

load("cpms_phenotype_merge.Rda")

ELO_genes2<-cpms_phenotype_merge%>%
  dplyr::select(genes,value,ELO_hrs,variable)%>%
  pivot_wider(names_from = "genes",values_from = "value")%>%
  dplyr::select(-variable)


#first, pick best gene which has the highest R2 on its own
TotalR2<-c()
  ELO_genes2$ELO_hrs_RANDOM<-ELO_genes2$ELO_hrs #THIS DOES NOT RANDOMIZE! Contrast with MultipleGenesPredictPheno_ELO.R
  ELO_genesRAN<-ELO_genes2[,-1]
  highestR2<-c()
  for (i in 1:7938){
    ELO_genes3<- ELO_genesRAN %>%
      dplyr::select(ELO_hrs_RANDOM, i)
    summary_fullModel3<-summary(lm(ELO_hrs_RANDOM~.,ELO_genes3))
    highestR2<-rbind(highestR2,c(colnames(ELO_genesRAN)[i],summary_fullModel3$r.squared))
  }
  
  highestR2<-as.data.frame(highestR2)
  colnames(highestR2)<-c("WBID","R2")
  highestR2$R2<-as.numeric(highestR2$R2)
  tail(highestR2[order(highestR2$R2),])
  tail(highestR2[order(highestR2$R2),])[6,1]
  R2_firstgene<-tail(highestR2[order(highestR2$R2),])[6,2]
  
  Gene<-c(tail(highestR2[order(highestR2$R2),])[6,1])
  for (j in 1:50){
    highestR2<-c()
    for (i in 1:7938){
      ELO_genes3<- ELO_genesRAN %>%
        dplyr::select(ELO_hrs_RANDOM, Gene[1:j],i)
      summary_fullModel3<-summary(lm(ELO_hrs_RANDOM~.,ELO_genes3))
      highestR2<-rbind(highestR2,c(colnames(ELO_genes3)[2+j],summary_fullModel3$r.squared))
      
    }
    highestR2<-as.data.frame(highestR2)
    highestR2$V2<-as.numeric(highestR2$V2)
    colnames(highestR2)<-c("Gene","FullR2")
    head(highestR2)
    highestR2 <- highestR2[order(-highestR2$FullR2),] 
    Gene<-c(Gene,highestR2[1,1])
    TotalR2<-rbind(TotalR2,c(0,j,Gene[1],R2_firstgene,highestR2[1,1],highestR2[1,2]))
  }



TotalR2<-as.data.frame(TotalR2)
write.table(TotalR2,"R2_shuffleddata_ELO_base.txt",quote = F,sep = "\t")


