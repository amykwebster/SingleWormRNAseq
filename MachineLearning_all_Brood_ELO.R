
library(glmnet)
load("cpms_phenotype_merge.Rda")
head(cpms_phenotype_merge)
unique(cpms_phenotype_merge$variable)

#BROOD - INCLUDING ALL 180 WORMS
cpms_phenotype_merge_brood<-cpms_phenotype_merge[,c(1,2,3,11)]
head(cpms_phenotype_merge_brood)

cpms_pheno_brood_wide<-pivot_wider(cpms_phenotype_merge_brood,names_from = "genes",values_from = "value",id_cols = "variable")
brood_sample_only<-subset(cpms_phenotype_merge,genes=="WBGene00000001")
head(brood_sample_only)
brood_sample_only<-brood_sample_only[,c(1,11)]
brood_sample_genes<-merge(brood_sample_only,cpms_pheno_brood_wide,by="variable")
head(brood_sample_genes[,1:4])

#dependent variable
dv<-brood_sample_genes$brood_size
pred<-as.matrix(brood_sample_genes[,3:7940])
head(pred[1:4,1:3])
pred<-scale(pred)

brood_predict_elastic<-c()
for (i in 1:180){
  set.seed(5)#added this in later..
  pred.train <- pred[-i,]
  dv.train <- dv[-i]
  pred.test <- pred[i,]
  dv.test <- dv[i]
  
  elastic<-cv.glmnet(pred.train, dv.train, type.measure="mse", alpha=0.5, family="gaussian")
  
  elastic.predicted <- predict(elastic, elastic$lambda.1se, new=pred.test)
  
  brood_predict_elastic<-rbind(brood_predict_elastic,c(brood_sample_genes[i,1],brood_sample_genes[i,2],dv.test,elastic.predicted[1]))
}
brood_predict_elastic

brood_predict_elastic<-as.data.frame(brood_predict_elastic)
brood_predict_elastic
colnames(brood_predict_elastic)<-c("sample","actual_brood1","actual_brood","predicted_brood")
class(brood_predict_elastic$actual_brood)
brood_predict_elastic$actual_brood<-as.numeric(brood_predict_elastic$actual_brood)
brood_predict_elastic$predicted_brood<-as.numeric(brood_predict_elastic$predicted_brood)

ggplot(brood_predict_elastic,aes(x=actual_brood,y=predicted_brood))+
  geom_point()+geom_smooth(method = "lm",formula = y~0+x)+
  ggtitle("Elastic Net ML, Leave-one-out, R2=0.53")+
  #geom_abline(slope=1,intercept = 0)+
  labs(x="Actual early brood",y="Predicted early brood")+
  theme_bw()+theme(aspect.ratio = 1)+xlim(70,160)+ylim(70,160)

summary(lm(predicted_brood~actual_brood,data=brood_predict_elastic)) #R2=0.53

ggplot(brood_predict_elastic,aes(x=diff))+
  geom_density()

head(brood_predict_elastic)
brood_predict_elastic$diff<-abs(brood_predict_elastic$predicted_brood - brood_predict_elastic$actual_brood)

write.table(brood_predict_elastic,"Brood_ML_Prediction_all.txt",quote = F,sep = "\t")



#ELO - INCLUDING ALL 180 WORMS
cpms_phenotype_merge_ELO<-cpms_phenotype_merge[,c(1,2,3,10)]
head(cpms_phenotype_merge_ELO)

cpms_pheno_ELO_wide<-pivot_wider(cpms_phenotype_merge_ELO,names_from = "genes",values_from = "value",id_cols = "variable")
ELO_sample_only<-subset(cpms_phenotype_merge,genes=="WBGene00000001")
head(ELO_sample_only)
ELO_sample_only<-ELO_sample_only[,c(1,10)]
ELO_sample_genes<-merge(ELO_sample_only,cpms_pheno_ELO_wide,by="variable")
head(ELO_sample_genes[,1:4])

#dependent variable
dv<-ELO_sample_genes$ELO_hrs
pred<-as.matrix(ELO_sample_genes[,3:7940])
head(pred[1:4,1:3])
pred<-scale(pred)

ELO_predict_elastic<-c()
for (i in 1:180){
  set.seed(5)#added this in later..
  pred.train <- pred[-i,]
  dv.train <- dv[-i]
  pred.test <- pred[i,]
  dv.test <- dv[i]
  
  elastic<-cv.glmnet(pred.train, dv.train, type.measure="mse", alpha=0.5, family="gaussian")
  
  elastic.predicted <- predict(elastic, elastic$lambda.1se, new=pred.test)
  
  ELO_predict_elastic<-rbind(ELO_predict_elastic,c(ELO_sample_genes[i,1],ELO_sample_genes[i,2],dv.test,elastic.predicted[1]))
}
ELO_predict_elastic

ELO_predict_elastic<-as.data.frame(ELO_predict_elastic)
ELO_predict_elastic
colnames(ELO_predict_elastic)<-c("sample","actual_ELO1","actual_ELO","predicted_ELO")
class(ELO_predict_elastic$actual_ELO)
ELO_predict_elastic$actual_ELO<-as.numeric(ELO_predict_elastic$actual_ELO)
ELO_predict_elastic$predicted_ELO<-as.numeric(ELO_predict_elastic$predicted_ELO)

summary(lm(predicted_ELO~actual_ELO,data=ELO_predict_elastic)) #R2=0.47
summary(lm(predicted_ELO~actual_ELO+0,data=ELO_predict_elastic)) #R2=0.9995


ggplot(ELO_predict_elastic,aes(x=actual_ELO,y=predicted_ELO))+
  geom_point()+
  geom_smooth(method = "lm")+
  #geom_smooth(method = "lm",formula = y~0+x)+
  ggtitle("Elastic Net ML, Leave-one-out, R2=0.48")+
  labs(x="Actual egg-laying onset",y="Predicted egg-laying onset")+
  theme_bw()+theme(aspect.ratio = 1)


head(ELO_predict_elastic)
ELO_predict_elastic$diff<-abs(ELO_predict_elastic$predicted_ELO - ELO_predict_elastic$actual_ELO)

write.table(ELO_predict_elastic,"ELO_ML_Prediction_all.txt",quote = F,sep = "\t")