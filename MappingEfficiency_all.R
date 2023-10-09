#Mapping efficiency for batches 1 and 2

#Each batch also includes 1 undetermined reads file, so 97 total
library(reshape2)
library(ggplot2)
library(edgeR)
library(nlme)


#BATCH 1
setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/singleworm_96samples/counts/")
#import all count files, then put them in a list. Assign original names to the objects in the list.
singlewormSummaryFiles= list.files(pattern= ".summary$",full.names = TRUE)
singlewormSummary=lapply(singlewormSummaryFiles,read.table, header=T)
melt(singlewormSummary[[1]])

singlewormMeltSummary<-list()
for (i in 1:97){
  singlewormMeltSummary[[i]]<-melt(singlewormSummary[[i]])
}
head(singlewormMeltSummary)

singlewormMeltDF<-do.call(rbind.data.frame,singlewormMeltSummary)

head(singlewormMeltDF,26)

singlewormTotalReads<-aggregate(singlewormMeltDF$value,list(library=singlewormMeltDF$variable),sum)
head(singlewormTotalReads)
singlewormTotalReads$library<-gsub("X20220222_singleworm.bams.","",as.character(singlewormTotalReads$library))
singlewormTotalReads$library<-gsub("_L001_R1_001.bam","",as.character(singlewormTotalReads$library))

setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/singleworm_96samples/")

head(singlewormTotalReads)
write.table(singlewormTotalReads,"TotalReads_singlewormlibraries.txt",quote = F,sep = "\t")

assignedReads<-subset(singlewormMeltDF,Status=="Assigned")
head(assignedReads)
assignedReads$variable<-gsub("X20220222_singleworm.bams.","",as.character(assignedReads$variable))
assignedReads$variable<-gsub("_L001_R1_001.bam","",as.character(assignedReads$variable))


singlewormTotalReads<-singlewormTotalReads[order(singlewormTotalReads$x),]
head(singlewormTotalReads)
singlewormTotalReads$library
ggplot(singlewormTotalReads,aes(x=library,y=x))+
  geom_bar(stat="identity")+theme_bw(base_size = 20)+labs(x="Library",y="Number of total reads")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_hline(yintercept = 4166667)+scale_x_discrete(limits=singlewormTotalReads$library)
singlewormTotalReads$library[c(1:95,97)]
#mapping efficiency
mappingEff<-merge(singlewormTotalReads,assignedReads,by.x = "library",by.y = "variable")
head(mappingEff)
tail(mappingEff)
mappingEff<-mappingEff[order(mappingEff$x),]

setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/")

head(mappingEff)
write.table(mappingEff,"MappingEfficiency_batch1.txt",quote = F,sep = "\t")

#BATCH 2

setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/singleworm_batch2/counts/")
#import all count files, then put them in a list. Assign original names to the objects in the list.
singlewormSummaryFiles_batch2= list.files(pattern= ".summary$",full.names = TRUE)
singlewormSummary_batch2=lapply(singlewormSummaryFiles_batch2,read.table, header=T)

singlewormMeltSummary_batch2<-list()
for (i in 1:97){
  singlewormMeltSummary_batch2[[i]]<-melt(singlewormSummary_batch2[[i]])
}

singlewormMeltDF<-do.call(rbind.data.frame,singlewormMeltSummary_batch2)
singlewormTotalReads<-aggregate(singlewormMeltDF$value,list(library=singlewormMeltDF$variable),sum)
singlewormTotalReads$library<-gsub("X20221220_singleworm_batch2.bams.","",as.character(singlewormTotalReads$library))
singlewormTotalReads$library<-gsub("_L002_R1_001.bam","",as.character(singlewormTotalReads$library))

setwd("/Users/amywebster/Documents/PhillipsLab/Aim1/singleworm_batch2/")

assignedReads<-subset(singlewormMeltDF,Status=="Assigned")
assignedReads$variable<-gsub("X20221220_singleworm_batch2.bams.","",as.character(assignedReads$variable))
assignedReads$variable<-gsub("_L002_R1_001.bam","",as.character(assignedReads$variable))

UnmappedReads<-subset(singlewormMeltDF,Status=="Unassigned_Unmapped")
UnmappedReads$variable<-gsub("X20221220_singleworm_batch2.bams.","",as.character(assignedReads$variable))
UnmappedReads$variable<-gsub("_L002_R1_001.bam","",as.character(assignedReads$variable))

singlewormTotalReads<-singlewormTotalReads[order(singlewormTotalReads$x),]

ggplot(singlewormTotalReads,aes(x=library,y=x))+
  geom_bar(stat="identity")+theme_bw(base_size = 10)+labs(x="Library",y="Number of total reads")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_hline(yintercept = 4166667)+
  geom_hline(yintercept = 1000000)+
  ggtitle("Total reads per library")+
  scale_x_discrete(limits=singlewormTotalReads$library)


mappingEff2<-merge(singlewormTotalReads,UnmappedReads,by.x = "library",by.y = "variable")
mappingEff2<-mappingEff2[order(mappingEff2$x),]

mappingEff2$PropMapped<-((mappingEff2$x - mappingEff2$value) / mappingEff2$x)

ggplot(mappingEff2,aes(x=library,y=PropMapped))+
  geom_bar(stat="identity")+theme_bw(base_size = 10)+labs(x="Library",y="Proportion of reads mapped")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_x_discrete(limits=singlewormTotalReads$library)+ggtitle("Mapping efficiency per library")


mappingEff<-merge(singlewormTotalReads,assignedReads,by.x = "library",by.y = "variable")
mappingEff<-mappingEff[order(mappingEff$x),]

mappingEff$PropMapped<-(mappingEff$value / mappingEff$x)

ggplot(mappingEff,aes(x=library,y=PropMapped))+
  geom_bar(stat="identity")+theme_bw(base_size = 10)+labs(x="Library",y="Proportion of reads assigned")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_x_discrete(limits=singlewormTotalReads$library)+
  ggtitle("Assigned reads")

write.table(mappingEff,"MappingEfficiency_batch2.txt",quote = F,sep = "\t")


#after exporting Batches 1 and 2 mapping effeciency, combine into csv file


