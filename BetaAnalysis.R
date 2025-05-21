setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/combined_analysis/")
library(gplots)
library(nlme)
library(lme4)
library(lmerTest)
load("cpms_phenotype_merge_8824.Rda")

#the following script will obtain Beta1 for all genes for both traits by looping through each gene and extracting 
    #the beta(s) from the model. 
#to obtain other Betas, alter fit_lmm and extract the relevant coefficients and normalize to Beta coefficients
#note that by looking at summary_fit_lmm$coefficients you can ensure that you extract the correct coefficients
#after extracting the coefficient, multiply by the standard deviation of the x variable and divide by the standard deviation of the y variable
#early brood Beta2, Beta3, and Beta4 can be extracted from this model:
#fit_lmm <- lmer(brood_size ~ value + ParentBinary + TempBinary + (1|replicate), Gene1_practice_AllBinary)
#ELO Beta2, Beta3, and Beta4 can be extracted from this model:
#fit_lmm <- lmer(ELO_hrs ~ value + ParentBinary + TempBinary + (1|replicate), Gene1_practice_AllBinary)
#Beta5 and Beta6 can be extracted from this model (not trait-dependent):
#fit_lmm <- lmer(value ~ ParentBinary + TempBinary + (1|replicate), Gene1_practice_AllBinary)

# Loop to get Beta1 for 8824 genes
gene_vector<-unique(cpms_phenotype_merge$genes)
df_mixedModel_Betas<-c()
for (i in 1:8824){
  new_row<-c()
  Gene1_practice<-subset(cpms_phenotype_merge,genes==gene_vector[i])
  Gene1_practice_YA<-subset(Gene1_practice,maternal_age=="YA")
  head(Gene1_practice_YA)
  Gene1_practice_OA<-subset(Gene1_practice,maternal_age=="OA")
  Gene1_practice_YA$ParentBinary<-"1"
  Gene1_practice_OA$ParentBinary<-"0"
  Gene1_practice_AllParentBinary<-rbind(Gene1_practice_YA,Gene1_practice_OA)
  Gene1_practice_constanttemp<-subset(Gene1_practice_AllParentBinary,environment=="20C")
  Gene1_practice_25C<-subset(Gene1_practice_AllParentBinary,environment=="25C_8hr")
  Gene1_practice_constanttemp$TempBinary<-"1"
  Gene1_practice_25C$TempBinary<-"0"
  Gene1_practice_AllBinary<-rbind(Gene1_practice_constanttemp,Gene1_practice_25C)
  fit_lmm <- lmer(brood_size ~ value + (1|replicate), Gene1_practice_AllBinary)
  summary_fit_lmm<-summary(fit_lmm)
  ExpBeta<-(summary_fit_lmm$coefficients[2,1]*sd(Gene1_practice_AllBinary$value))/sd(Gene1_practice_AllBinary$brood_size)
  ExpPval<-summary_fit_lmm$coefficients[2,5]
  new_row<-c(gene_vector[i],ExpBeta,ExpPval)
  df_mixedModel_Betas<-rbind(df_mixedModel_Betas,new_row)
}

df_mixedModel_Betas<-as.data.frame(df_mixedModel_Betas)
head(df_mixedModel_Betas)
colnames(df_mixedModel_Betas)<-c("gene","exp_beta","exp_pval")
df_mixedModel_Betas<-df_mixedModel_Betas%>%
  mutate_at(2:3,as.numeric)

as_tibble(df_mixedModel_Betas)
range(abs(df_mixedModel_Betas$exp_beta))




#Loop to get Beta1 for ELO for 8824 genes
gene_vector<-unique(cpms_phenotype_merge$genes)
df_mixedModel_Betas2<-c()
for (i in 1:8824){
  new_row<-c()
  Gene1_practice<-subset(cpms_phenotype_merge,genes==gene_vector[i])
  #Gene1_practice<-subset(cpms_phenotype_merge,genes=="WBGene00019611")
  Gene1_practice_YA<-subset(Gene1_practice,maternal_age=="YA")
  head(Gene1_practice_YA)
  Gene1_practice_OA<-subset(Gene1_practice,maternal_age=="OA")
  Gene1_practice_YA$ParentBinary<-"1"
  Gene1_practice_OA$ParentBinary<-"0"
  Gene1_practice_AllParentBinary<-rbind(Gene1_practice_YA,Gene1_practice_OA)
  Gene1_practice_constanttemp<-subset(Gene1_practice_AllParentBinary,environment=="20C")
  Gene1_practice_25C<-subset(Gene1_practice_AllParentBinary,environment=="25C_8hr")
  Gene1_practice_constanttemp$TempBinary<-"1"
  Gene1_practice_25C$TempBinary<-"0"
  Gene1_practice_AllBinary<-rbind(Gene1_practice_constanttemp,Gene1_practice_25C)
  fit_lmm <- lmer(ELO_hrs ~ value + (1|replicate), Gene1_practice_AllBinary)
  summary_fit_lmm<-summary(fit_lmm)
  ExpBeta<-(summary_fit_lmm$coefficients[2,1]*sd(Gene1_practice_AllBinary$value))/sd(Gene1_practice_AllBinary$ELO_hrs)
  ExpPval<-summary_fit_lmm$coefficients[2,5]
  new_row<-c(gene_vector[i],ExpBeta,ExpPval)
  df_mixedModel_Betas2<-rbind(df_mixedModel_Betas2,new_row)
}

df_mixedModel_Betas2<-as.data.frame(df_mixedModel_Betas2)
head(df_mixedModel_Betas2)
colnames(df_mixedModel_Betas2)<-c("gene","exp_beta","exp_pval")
df_mixedModel_Betas2<-df_mixedModel_Betas2%>%
  mutate_at(2:3,as.numeric)

as_tibble(df_mixedModel_Betas2)
range(abs(df_mixedModel_Betas2$exp_beta))


