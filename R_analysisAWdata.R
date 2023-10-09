
#Script containing analysis by Peter Sarkies 

load("Ahringer_ref_table.Rdata")

varG<-read.table("Zscores_MixedModel2.csv", sep=",", header=T, stringsAsFactors=F)

Act<-rep(0, length=nrow(varG))
Reg<-rep(0, length=nrow(varG))
Germline<-rep(0, length=nrow(varG))
Somatic<-rep(0, length=nrow(varG))
Ubiq<-rep(0, length=nrow(varG))

ActGenes<-Ahringer_ref_table$Gene_mapped[which(Ahringer_ref_table$Chromatin_domain=="A")]
RegGenes<-Ahringer_ref_table$Gene_mapped[which(Ahringer_ref_table$Chromatin_domain=="R")]

GermlineGenes<-Ahringer_ref_table$Gene_mapped[grep("Germline", Ahringer_ref_table$Tissue_specificity)]
SomaticClasses<-c("Neurons","Muscle","Intest","Hypod","Soma")
SomaticGenes<-c()
for(i in 1:5){

SomaticGenes<-c(SomaticGenes,Ahringer_ref_table$Gene_mapped[grep(SomaticClasses[i],Ahringer_ref_table$Tissue_specificity)])
}


UbiqGenes<-Ahringer_ref_table$Gene_mapped[grep("Ubiq",Ahringer_ref_table$Tissue_specificity)]

Act[which(varG$symbol%in%ActGenes==T|varG$sequence%in%ActGenes==T)]<-1
Reg[which(varG$symbol%in%RegGenes==T|varG$sequence%in%RegGenes==T)]<-1

Germline[which(varG$symbol%in%GermlineGenes==T|varG$sequence%in%GermlineGenes==T)]<-1

Somatic[which(varG$symbol%in%SomaticGenes==T|varG$sequence%in%SomaticGenes==T)]<-1
Ubiq[which(varG$symbol%in%UbiqGenes==T|varG$sequence%in%UbiqGenes==T)]<-1

varZGlm<-glm((formula = varG$zscore_controlled_var ~ Act + Reg + Germline + Ubiq + 
    Somatic))
    
PlotData<-barplot(varZGlm$coefficients[-1],ylim=c(-0.1,0.5), ylab="effect (combinedZ)")

pVals<-summary(varZGlm)$coefficients[-1,4]
pVals<-format(pVals, scientific=T,digits=3)

text(PlotData[,1],varZGlm$coefficients[-1]+0.09,paste("p=",pVals,sep=""))
dev.copy(pdf, "VarGenes_glmCoeff_ud.pdf")
dev.off()


boxplot(varG$zscore_controlled_var[which(Reg==1)],varG$zscore_controlled_var[which(Act==1)],notch=T,names=c("Regulated\nH3K27me3","Active\nH3K36me3"),col=c("red","blue"),outline=F,ylab="Zscore_combined")

wilcox.test(varG$zscore_controlled_var[which(Reg==1)],varG$zscore_controlled_var[which(Act==1)])$p.value

text(1.4,2,"p=5.94e-52")
dev.copy(pdf, "RegvAct_combinedZ_comparison_ud.pdf")
dev.off()


glmList<-list()

for(i in 4:8){glmList[[i-3]]<-glm(formula=varG[,i]~Act+Reg+Germline+Ubiq+Somatic)}

colPlot<-c("blue","red","orange","grey","green")
plot(1:5,glmList[[1]]$coefficients[-1],pch=16,col=colPlot,ylim=c(-0.5,0.5),ylab="Effect on Z score",xaxt="null")
axis(side=1, at=c(1,2,3,4,5),c("Act","Reg","Germline","Ubiq.","Soma"))
for(i in 2:5){points(1:5,glmList[[i]]$coefficients[-1],pch=i-1,col=colPlot ,ylim=c(-0.5,0.5))}
legend("bottomright", colnames(varG)[4:8],pch=c(16,1,2,3,4),col="black")

brood<-read.table("BroodGenes_avgCPMS.csv", sep=",", header=T, stringsAsFactors=F)

ActB<-rep(0, length=nrow(brood))
RegB<-rep(0, length=nrow(brood))
GermlineB<-rep(0, length=nrow(brood))
SomaticB<-rep(0, length=nrow(brood))
UbiqB<-rep(0, length=nrow(brood))


ActB[which(brood$symbol%in%ActGenes==T|brood$sequence%in%ActGenes==T)]<-1
RegB[which(brood$symbol%in%RegGenes==T|brood$sequence%in%RegGenes==T)]<-1

GermlineB[which(brood$symbol%in%GermlineGenes==T|brood$sequence%in%GermlineGenes==T)]<-1

SomaticB[which(brood$symbol%in%SomaticGenes==T|brood$sequence%in%SomaticGenes==T)]<-1
UbiqB[which(brood$symbol%in%UbiqGenes==T|brood$sequence%in%UbiqGenes==T)]<-1


BroodSOMA<-brood$CoeffCPMNorm[which(SomaticB==1)]
RegSOMA<-RegB[which(SomaticB==1)]
ActSOMA<-ActB[which(SomaticB==1)]

broodCPMglmSoma<-glm(BroodSOMA~RegSOMA+ActSOMA)

BroodGERM<-brood$CoeffCPMNorm[which(GermlineB==1)]
RegGERM<-RegB[which(GermlineB==1)]
ActGERM<-ActB[which(GermlineB==1)]

broodCPMglmGerm<-glm(BroodGERM~RegGERM+ActGERM)



broodCPMglm<-glm(brood$CoeffCPMNorm~ActB+RegB+GermlineB+UbiqB+SomaticB)
PlotData2<-barplot(broodCPMglm$coefficients[-1],names=c("Active","Regulated","Germline","Ubiquitous","Somatic"),ylab="broodSize_effect",ylim=c(-17,17))

pVals2<-summary(broodCPMglm)$coefficients[-1,4]
pVals2<-format(pVals2,scientific=T,digits=3)
PlotPoints<-broodCPMglm$coefficients[-1]
PlotPoints[PlotPoints<0]<-PlotPoints[PlotPoints<0]-0.5
PlotPoints[PlotPoints>0]<-PlotPoints[PlotPoints>0]+0.5
text(PlotData[,1],PlotPoints,paste("p=",pVals2,sep=""))

dev.copy(pdf, "effect_size_all_brood.pdf")
dev.off()



Glm_resultsALL<-c(broodCPMglm$coefficients[-1],broodCPMglmSoma$coefficients[-1],broodCPMglmGerm$coefficients[-1])
PvalALL<-c(summary(broodCPMglm)$coefficients[-1,4],summary(broodCPMglmSoma)$coefficients[-1,4],summary(broodCPMglmGerm)$coefficients[-1,4])

CombPlot<-barplot(Glm_resultsALL, las=2, names=c("Active all\nH3K36me3","Regulated all\nH3K27me3","Germline all","Ubiq. all","Soma all","Regulated\nsoma","Act\nsoma","Regulated\ngermline","Active\ngermline"),ylab="effect on brood",ylim=c(-16,15))

PlotPointsComb<-Glm_resultsALL
PlotPointsComb[which(PlotPointsComb>0)]<-PlotPointsComb[which(PlotPointsComb>0)]+0.35
PlotPointsComb[which(PlotPointsComb<0)]<-PlotPointsComb[which(PlotPointsComb<0)]-0.45
PvalALL<-format(PvalALL, scientific=T, digits=2)
PvalALL<-p.adjust(PvalALL,method="bonferroni")
text(CombPlot,PlotPointsComb,paste("p=",PvalALL,sep=""),cex=0.5)







