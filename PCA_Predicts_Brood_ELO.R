#Analysis to determine how well gene expression predicts traits
library(reshape2)
library(tidyverse)
library(edgeR)
library(gridExtra)
setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/combined_analysis/")


setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/combined_analysis/counts/")
singlewormCountFiles= list.files(pattern= ".txt$",full.names = TRUE)
singlewormCounts=lapply(singlewormCountFiles,read.table, header=T)
setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/singleworm_96samples/")
WS273_geneNames<-read.csv("WS273_geneNames.csv",header = T) #this has the extra information beyond 

setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/combined_analysis/")
singleworm_phenotypes<-read.csv("singleworm_96samples_combined.csv",header = T)
singleworm_phenotypes<-singleworm_phenotypes[,1:11]

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
dim(d2)[1] #8824 genes, 7938 for the >1 in all 180 filter
cpm_d2<-cpm(d2,normalized.lib.sizes = TRUE) #make a counts per million object containing normalized CPM
cpm_d2<-as.data.frame(cpm_d2)
head(cpm_d2)

head(cpm_d2)
cpm_d2_df<-data.frame(cpm_d2)
cpm_d2_df$mean<-rowMeans(cpm_d2_df) #calculate mean CPM across all libraries
#change the 1:10 depending on how many libraries there are
cpm_d2_df2<-cpm_d2_df[,1:180]/cpm_d2_df$mean #mean normalize
cpm_d2_df2<-log2(cpm_d2_df2+1) #log2 transform 
pca = prcomp(t(cpm_d2_df2)) #principal component analysis (PCA) on the log2 mean normalized CPM values
summary_pca<-summary(pca)
summary_pca

pca_genes<-pca$x
pca_genes_dataframe<-as.data.frame(pca_genes)
head(pca_genes_dataframe)
pca_genes_dataframe$sample<-rownames(pca_genes_dataframe)
library(stringr)
pca_genes_dataframe$sample <- str_extract(pca_genes_dataframe$sample, "(?<=X)\\d+(?=\\_S)")

pca_genes_dataframe2<-merge(pca_genes_dataframe,singleworm_phenotypes,by.x = "sample",by.y = "RNA_sample_number")
head(pca_genes_dataframe2)
pca_genes_dataframe2<-pca_genes_dataframe2[,c(-1,-182:-188,-190:-191)]
pca_genes_dataframe3<-merge(pca_genes_dataframe,singleworm_phenotypes,by.x = "sample",by.y = "RNA_sample_number")
pca_genes_dataframe3_ELO<-pca_genes_dataframe3[,c(-1,-182:-187,-189:-191)]
head(pca_genes_dataframe3_ELO)


random_brood<-c()
for (i in 1:100){
  set.seed(i)
  rand_order<-sample(nrow(pca_genes_dataframe2))
  pca_genes_dataframe2$brood_RANDOM<-pca_genes_dataframe2$brood_size[rand_order]
  for (j in 2:178){
    summary_random<-summary(lm(brood_RANDOM ~ ., data = pca_genes_dataframe2[,c(1:j,182)]))
    summary_random$r.squared
    random_brood<-rbind(random_brood,c(i,j,summary_random$r.squared))
  }
}

for (j in 2:178){
  summary_random<- summary(lm(brood_size ~ ., data = pca_genes_dataframe2[,c(1:j,181)]))
  summary_random$r.squared
  random_brood<-rbind(random_brood,c(0,j,summary_random$r.squared))
}
tail(random_brood)

random_brood<-as.data.frame(random_brood)
colnames(random_brood)<-c("Trial","PCs","R2")
head(random_brood)
class(random_brood$PCs)
class(random_brood$R2)

head(random_brood)
unique(random_brood$Trial)

random_brood_only<-random_brood%>%
  filter(Trial!="0")%>%
  group_by(PCs)%>%
  summarize(Average=mean(R2),StDev=sd(R2))

random_brood_only

EarlyBrood_PCs<-ggplot(random_brood_only,aes(x=PCs,y=Average))+
  geom_line(alpha=0.4,size=0.7)+
  #xlim(0,100)+
  geom_ribbon(aes(ymin=Average-StDev,ymax=Average+StDev),alpha=0.2)+
  geom_line(data=subset(random_brood,Trial=="0"),aes(x=PCs,y=R2),color="red")+
  theme(aspect.ratio = 1)+ggtitle("Prediction of early brood by PCs")+
  theme_bw(base_size = 12)+theme(aspect.ratio = 1)+
  labs(x="PCs included in multiple regression",y="Total R2 of multiple regression")

EarlyBrood_PCs

random_ELO<-c()
for (i in 1:100){
  set.seed(i)
  rand_order<-sample(nrow(pca_genes_dataframe3_ELO))
  pca_genes_dataframe3_ELO$ELO_RANDOM<-pca_genes_dataframe3_ELO$ELO_hrs[rand_order]
  for (j in 2:178){
    summary_random<-summary(lm(ELO_RANDOM ~ ., data = pca_genes_dataframe3_ELO[,c(1:j,182)]))
    summary_random$r.squared
    random_ELO<-rbind(random_ELO,c(i,j,summary_random$r.squared))
  }
}

for (j in 2:178){
  summary_random<- summary(lm(ELO_hrs ~ ., data = pca_genes_dataframe3_ELO[,c(1:j,181)]))
  summary_random$r.squared
  random_ELO<-rbind(random_ELO,c(0,j,summary_random$r.squared))
}
tail(random_ELO)

random_ELO<-as.data.frame(random_ELO)
colnames(random_ELO)<-c("Trial","PCs","R2")
head(random_ELO)
class(random_ELO$PCs)
class(random_ELO$R2)

head(random_ELO)

random_ELO_only<-random_ELO%>%
  filter(Trial!="0")%>%
  group_by(PCs)%>%
  summarize(Average=mean(R2),StDev=sd(R2))

random_ELO_only

ELO_PCs<-ggplot(random_ELO_only,aes(x=PCs,y=Average))+
  geom_line(alpha=0.4,size=0.7)+
  #xlim(0,100)+
  geom_ribbon(aes(ymin=Average-StDev,ymax=Average+StDev),alpha=0.2)+
  geom_line(data=subset(random_ELO,Trial=="0"),aes(x=PCs,y=R2),color="red")+
  theme(aspect.ratio = 1)+ggtitle("Prediction of ELO by PCs")+
  theme_bw(base_size = 12)+theme(aspect.ratio = 1)+
  labs(x="PCs included in multiple regression",y="Total R2 of multiple regression")


subset(random_ELO,Trial=="0")


grid.arrange(EarlyBrood_PCs,ELO_PCs,ncol=2)


