#sample script for processing coverage files: singlewormCoverageFiles=list.files(pattern="") is changed to process the entire genome for mutation calculations
#set working directory & load libraries
setwd("~/phillipslab/singleworm_mutationAnalysis/bams_all")
library(tidyverse)
#BiocManager::install("GenomicAlignments")
library(GenomicAlignments)

#bring in bam file and coverage file
singlewormCoverageFiles= list.files(pattern= "_chrI.coverage$",full.names = TRUE)
singlewormBams=list.files(pattern = ".bam$",full.names = TRUE)

NumberNucleotides<-c()
Unambiguous_ALL<-c()
Ambiguous_ALL<-c()
MultiGeno_All<-c()

#REPLACE with 1:5 to start..instead of length(singlewormCoverageFiles)
for (i in 1:length(singlewormCoverageFiles)){
  coverageFile<-read.table(singlewormCoverageFiles[i],header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
  colnames(coverageFile)<-c("chr","coord","coverage")
  #coverageFile_20reads<-subset(coverageFile,coverage>20 coord<20000)
  #coverageFile_20reads<-subset(coverageFile,coverage>20 & coord>19999 & coord<1000000)
  #coverageFile_20reads<-subset(coverageFile,coverage>20 & coord>999999 & coord<3000000)
  #coverageFile_20reads<-subset(coverageFile,coverage>20 & coord>2999999 & coord<4000000)
  #coverageFile_20reads<-subset(coverageFile,coverage>20 & coord>3999999 & coord<6000000)
  #coverageFile_20reads<-subset(coverageFile,coverage>20 & coord>5999999 & coord<8000000)
  #coverageFile_20reads<-subset(coverageFile,coverage>20 & coord>7999999 & coord<10000000)
  #coverageFile_20reads<-subset(coverageFile,coverage>20 & coord>9999999 & coord<12000000)
  #coverageFile_20reads<-subset(coverageFile,coverage>20 & coord>11999999 & coord<14000000)
  coverageFile_20reads<-subset(coverageFile,coverage>20 & coord>13999999)
  
  if ( dim(coverageFile_20reads)[1] > 0){
  bamFile1<-BamFile(singlewormBams[i])
  coverageFile_20reads$coordForNucPiles<-paste0(coverageFile_20reads$chr,":",coverageFile_20reads$coord,"-",coverageFile_20reads$coord)
  coverageFile_20reads_Nucleotides<-coverageFile_20reads$coordForNucPiles
  my_GPOI <- GPos((coverageFile_20reads_Nucleotides)[1:length(coverageFile_20reads_Nucleotides)])
  seqinfo(my_GPOI) <- merge(seqinfo(my_GPOI), seqinfo(bamFile1))
  seqlevels(my_GPOI) <- seqlevelsInUse(my_GPOI)
  
  ## Load the BAM file in a GAlignments object. Note that we load only
  ## the reads aligned to the sequences in 'seqlevels(my_GPOI)'. Also,
  ## in order to be consistent with applyPileups() and SAMtools (m)pileup,
  ## we filter out the following BAM records:
  ##   - secondary alignments (flag bit 0x100);
  ##   - reads not passing quality controls (flag bit 0x200);
  ##   - PCR or optical duplicates (flag bit 0x400).
  ## See ?ScanBamParam and the SAM Spec for more information. 
  
  which <- as(seqinfo(my_GPOI), "GRanges")
  flag <- scanBamFlag(isSecondaryAlignment=FALSE,
                      isNotPassingQualityControls=FALSE,
                      isDuplicate=FALSE)
  what <- c("seq", "qual")
  param <- ScanBamParam(flag=flag, what=c("seq", "qual"), which=which)
  gal <- readGAlignments(bamFile1, param=param)
  seqlevels(gal) <- seqlevels(my_GPOI) 
  
  scanbam_param <- ScanBamParam(flag=flag, which=my_GPOI)
  pileup_param <- PileupParam(max_depth=5000,
                              min_base_quality=0,
                              distinguish_strands=FALSE)
  piled_nucs<-pileup(bamFile1, scanBamParam=scanbam_param, pileupParam=pileup_param)
  
  
  UnambiguousGeno<-piled_nucs %>%
    group_by(seqnames,pos,which_label)%>%
    filter(n()==1)%>%
    mutate(genotype_class="unambiguous",sampleBam=singlewormBams[i],sampleCoverage=singlewormCoverageFiles[i])
  
  Unambiguous_ALL<-rbind(Unambiguous_ALL,UnambiguousGeno)
  
  AmbiguousGeno<-piled_nucs %>%
    group_by(seqnames,pos,which_label)%>%
    filter(n()==2)%>%
    mutate(genotype_class="ambiguous_2geno",sampleBam=singlewormBams[i],sampleCoverage=singlewormCoverageFiles[i])
  Ambiguous_ALL<-rbind(Ambiguous_ALL,AmbiguousGeno)
  
  MoreThanTwoGeno<-piled_nucs %>%
    group_by(seqnames,pos,which_label)%>%
    filter(n()>2)%>%
    mutate(genotype_class="ambiguous_morethan2",sampleBam=singlewormBams[i],sampleCoverage=singlewormCoverageFiles[i])
  MultiGeno_All<-rbind(MultiGeno_All,MoreThanTwoGeno)
  
  NumberNucleotides<-rbind(NumberNucleotides,c(singlewormCoverageFiles[i],dim(coverageFile_20reads)[1],median(coverageFile$coverage),dim(UnambiguousGeno)[1],dim(AmbiguousGeno)[1],dim(MoreThanTwoGeno)[1]))
  
}

}

colnames(NumberNucleotides)<-c("CoverageFileName","NucsWith20Reads","MedianCoverage","UnambiguousNucs","AmbiguousNucs","MoreThanTwoGenos")

PropGeno<-Unambiguous_ALL%>%
  group_by(seqnames,pos,which_label,nucleotide)%>%
  summarize(WormsWithNucleotide=n())%>% #211 #170 have more than 1, 41 have 1 (170x2 + 41 makes sense) in my pilot
  group_by(which_label)%>%
  mutate(WormsWithSite=sum(WormsWithNucleotide))%>%
  mutate(PropWithSameNuc=WormsWithNucleotide/WormsWithSite)%>%
  group_by(seqnames,pos,which_label,nucleotide,WormsWithSite)%>%
  summarize(MaxProp=max(PropWithSameNuc))

PropGeno_candidateMuts<-PropGeno%>%
  filter(MaxProp<1)


write.table(PropGeno,"GenotypeProportions_Unambiguous_chrI.txt",quote = F,sep = "\t")
write.table(PropGeno_candidateMuts,"GenotypeProportions_UnambiguousMuts_chrI.txt",quote = F,sep = "\t")

write.table(Unambiguous_ALL,"Unambiguous_ALL_data_chrI.txt",quote = F,sep = "\t")
write.table(Ambiguous_ALL,"Ambiguous_ALL_data_chrI.txt",quote = F,sep = "\t")
write.table(MultiGeno_All,"MultiGeno_ALL_data_chrI.txt",quote = F,sep = "\t")

write.table(NumberNucleotides,"NucleotideSummaryStats_chrI.txt",quote = F,sep = "\t")
