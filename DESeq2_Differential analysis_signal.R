
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(GenomicRanges)
library(chromVAR) ## For FRiP analysis and differential analysis
library(DESeq2) ## For differential analysis section
#library(ggpubr) ## For customizing figures
#library(corrplot) ## For correlation plot

# Run one clone at the time!

# Use a tsv files for samples annotation

samples_anno <- read.table("Clone1.txt", header = TRUE, sep = '\t')
print(samples_anno)
 
samples = unique(samples_anno$Sample)
reps = unique(samples_anno$Rep)

mPeak = GRanges()
for(samp in samples){
  sample_specific = samples_anno %>% filter(Sample == samp)
  for(rep in sample_specific$Rep){
    sprintf("Running for sample %s rep %s", samp, rep)
    rep_specific = sample_specific %>% filter(Rep == rep)
    path = paste0("/path/SEACR/", rep_specific$Filename, "_seacr_top0.05.peaks.stringent.bed")
    peakRes = read.table(path, header = FALSE, fill = TRUE)
    mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)}
} 
# Create a master peak list merging all the peaks called for each sample
masterPeak = reduce(mPeak)
#write.table(masterPeak, file="/cluster/projects/lupiengroup/People/ggrillo/Cut_and_Run/DESeq2/Differential_signal_intensity_H3K27ac_CHECK1.txt", sep="\t", quote=F, row.names=F)

# Get the fragment counts for each peak in the master peak list

#bamDir = "/cluster/projects/lupiengroup/People/ankita/giacomo/cutrun/map3/alignment/bam"
countMat = matrix(NA, length(masterPeak), length(samples)*length(reps))
## overlap with bam file to get count
i = 1
for(samp in samples){
  sample_specific = samples_anno %>% filter(Sample == samp)
  for(rep in sample_specific$Rep){
    sprintf("Running for sample %s rep %s", samp, rep)
    rep_specific = sample_specific %>% filter(Rep == rep)
    bamFile = paste0("/path/alignment/bam/", rep_specific$Filename, "_bowtie2.mapped.bam")
    #bamFile = paste0(bamDir, "/", hist, "_bowtie2.mapped.bam")
    fragment_counts <- getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
    countMat[, i] = counts(fragment_counts)[,1]
    i = i + 1
  }
}
colnames(countMat) = paste(rep(samples, 2), rep(reps, each = 2), sep = "_")
#write.table(countMat, file="/cluster/projects/lupiengroup/People/ggrillo/Cut_and_Run/DESeq2/Differential_signal_intensity_H3K27ac_CHECK2.txt", sep="\t", quote=F, row.names=F)
# Sequencing depth normalization and differential enriched peaks detection
rownames(countMat)=paste(as.character(seqnames(masterPeak)),as.character(start(masterPeak)),as.character(end(masterPeak)),sep="_")
selectR = which(rowSums(countMat) > 5) ## remove low count genes
#print(selectR)
dataS = countMat[selectR,]
#print(dataS)
print("rownames")
print(head(rownames(dataS)))
r<-rownames(dataS)

condition = factor(rep(samples, each = length(reps)))
dds = DESeqDataSetFromMatrix(countData = dataS,
                              colData = DataFrame(condition),
                              design = ~ condition)
DDS = DESeq(dds)
normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")
res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")
print(res)

countMatDiff = cbind(data.frame(dataS), data.frame(normDDS), data.frame(res))
rownames(countMatDiff)<-r
write.table(countMatDiff, file="Differential_signal_intensity.txt", sep="\t", quote=F, row.names=T)

