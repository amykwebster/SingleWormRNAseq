library(lme4)

load("cpms_phenotype_merge.Rda")

gene_vector<-unique(cpms_phenotype_merge$genes)

DE_residual_analysis<-c()
for (i in 4001:7938){
  tryCatch({
    Gene1_practice<-subset(cpms_phenotype_merge,genes==gene_vector[i])
    Gene1_model<-glmer.nb(value~maternal_age+environment+(1|replicate),data=Gene1_practice)
    Gene1Summary<-summary(Gene1_model)
    Gene1Residuals<-as.data.frame(Gene1Summary$residuals)
    Gene1Residuals$SquaredRes<-Gene1Summary$residuals*Gene1Summary$residuals
    DE_residual_analysis<-rbind(DE_residual_analysis,c(gene_vector[i],mean(Gene1_practice$value),Gene1Summary$coefficients[2,4],Gene1Summary$coefficients[3,4],sum(Gene1Residuals$SquaredRes),sum((fitted(Gene1_model)-Gene1_practice$value)^2), sum((mean(Gene1_practice$value)-Gene1_practice$value)^2)))
  },error=function(e){})
}

DE_residual_analysis<-as.data.frame(DE_residual_analysis)
DE_residual_analysis
colnames(DE_residual_analysis)<-c("genes","mean_CPM","p_parentage","p_temp","resid_ss","calc_resid_ss","calc_total_ss")
DE_residual_analysis$mean_CPM<-as.numeric(DE_residual_analysis$mean_CPM)
DE_residual_analysis$p_parentage<-as.numeric(DE_residual_analysis$p_parentage)
DE_residual_analysis$p_temp<-as.numeric(DE_residual_analysis$p_temp)
DE_residual_analysis$resid_ss<-as.numeric(DE_residual_analysis$resid_ss)
DE_residual_analysis$calc_resid_ss<-as.numeric(DE_residual_analysis$calc_resid_ss)
DE_residual_analysis$calc_total_ss<-as.numeric(DE_residual_analysis$calc_total_ss)
DE_residual_analysis$unexplained_var<-(DE_residual_analysis$calc_resid_ss / DE_residual_analysis$calc_total_ss)


write.table(DE_residual_analysis,"DiffExp_ResidAnalysis_SecondHalf.txt",quote = F,sep = "\t")
