#!/bin/bash

histName=$1 #"K27me3"

module load bowtie2/2.4.5
module load samtools/1.14
module load bedtools/2.27.1
module load fastqc/0.11.5
module load MACS/2.2.5

cores=5
projPath="/cluster/projects/path/"
fastqdir="/cluster/projects/fastqpath/"
seacr="${projPath}/SEACR/SEACR_1.3.sh"
trimmomaticpath="${projPath}/trimmomatic-0.36.jar"
ref="/cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/Bowtie2Index/genome"
spikeInRef="/cluster/tools/data/genomes/yeast/sacCer3/iGenomes/Sequence/Bowtie2Index/genome"
chromSize="${projPath}/hg38.chrom.sizes"

# ## 0. merge technical replicate
# echo "0. Merge the technical replicate"
# mkdir -p ${projPath}/fastq
# cat ${projPath}/data/${histName}/*_R1_*.fastq.gz >${projPath}/fastq/${histName}_R1.fastq.gz
# cat ${projPath}/data/${histName}/*_R2_*.fastq.gz >${projPath}/fastq/${histName}_R2.fastq.gz


## 0. Trim
echo "1. Trim"
mkdir -p ${projPath}/trim/${histName}
java -jar ${trimmomaticpath} PE -threads $cores ${fastqdir}/${histName}_R1_001.fastq.gz ${fastqdir}/${histName}_R2_001.fastq.gz -baseout ${projPath}/trim/${histName}/${histName}.fastq.gz ILLUMINACLIP:/cluster/projects/lupiengroup/pipelines/CUT-RUNTools2/CUT-RUNTools-2.0/adapters/Truseq3.PE.fa:2:15:4:1:true MINLEN:20

## 1. fastqc
echo "0. FastQC"
mkdir -p ${projPath}/fastqFileQC/${histName}

fastqc -o ${projPath}/fastqFileQC/${histName} -f fastq ${projPath}/trim/${histName}/${histName}_1P.fastq.gz
fastqc -o ${projPath}/fastqFileQC/${histName} -f fastq ${projPath}/trim/${histName}/${histName}_2P.fastq.gz


## 2. bowtie2 alignment
echo "2. Bowtie2 Alignment"


mkdir -p ${projPath}/alignment/sam/bowtie2_summary
mkdir -p ${projPath}/alignment/bam
mkdir -p ${projPath}/alignment/bed
mkdir -p ${projPath}/alignment/bedgraph

# ## bowtie2-build path/to/hg38/fasta/hg38.fa /path/to/bowtie2Index/hg38
bowtie2 -q -I 10 -X 700  --local --very-sensitive-local --no-mixed --no-discordant --no-unal --phred33 -p ${cores} -x ${ref} -1 ${projPath}/trim/${histName}/${histName}_1P.fastq.gz -2 ${projPath}/trim/${histName}/${histName}_2P.fastq.gz -S ${projPath}/alignment/sam/${histName}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt

mkdir -p $projPath/alignment/sam/fragmentLen

# ## Extract the 9th column from the alignment sam file which is the fragment length
samtools view -F 0x04 $projPath/alignment/sam/${histName}_bowtie2.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >$projPath/alignment/sam/fragmentLen/${histName}_fragmentLen.txt

## 3. In preparation for SEACR peak calling
echo "3. Data cleaning and format convertion"
samtools view -bS -F 0x04 $projPath/alignment/sam/${histName}_bowtie2.sam >$projPath/alignment/bam/${histName}_bowtie2.mapped.bam

bedtools bamtobed -i $projPath/alignment/bam/${histName}_bowtie2.mapped.bam -bedpe >$projPath/alignment/bed/${histName}_bowtie2.bed
awk '$1==$4 && $6-$2 < 1000 {print $0}' $projPath/alignment/bed/${histName}_bowtie2.bed >$projPath/alignment/bed/${histName}_bowtie2.clean.bed
cut -f 1,2,6 $projPath/alignment/bed/${histName}_bowtie2.clean.bed | grep '^chr[0-9XY]*'$'\t'  | sort -k1,1 -k2,2n -k3,3n  >$projPath/alignment/bed/${histName}_bowtie2.fragments.bed

## replicate reproducibility analysis
#binLen=500
#awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $projPath/alignment/bed/${histName}_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >$projPath/alignment/bed/${histName}_bowtie2.fragmentsCount.bin$binLen.bed


## 4. Spike-in calibration
echo "4. Spike-in calibration"


# ## bowtie2-build path/to/Ecoli/fasta/Ecoli.fa /path/to/bowtie2Index/Ecoli
bowtie2 -q -I 10 -X 700  --local --very-sensitive-local --no-mixed --no-discordant --no-overlap --no-dovetail --no-unal --phred33 -p ${cores} -x ${spikeInRef} -1 ${projPath}/trim/${histName}/${histName}_1P.fastq.gz -2 ${projPath}/trim/${histName}/${histName}_2P.fastq.gz -S $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam &> $projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.txt

seqDepthDouble=`samtools view -F 0x04 $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam | wc -l`
seqDepth=$((seqDepthDouble/2))
echo $seqDepth >$projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.seqDepth


if [[ "$seqDepth" -gt "1" ]]; then
     scale_factor=`echo "10000 / $seqDepth" | bc -l`
     echo "Scaling factor for $histName is: $scale_factor!"
     bedtools genomecov -bg -scale $scale_factor -i $projPath/alignment/bed/${histName}_bowtie2.fragments.bed -g $chromSize >$projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph
fi


## 5. SEACR peak calling

# histControl=$2
mkdir -p $projPath/peakCalling/SEACR

# bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
#      $projPath/alignment/bedgraph/${histControl}_bowtie2.fragments.normalized.bedgraph \
#      non stringent $projPath/peakCalling/SEACR/${histName}_seacr_control.peaks

bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph 0.05 non stringent $projPath/peakCalling/SEACR/${histName}_seacr_top0.05.peaks

## 6. macs2 peak calling
mkdir -p $projPath/peakCalling/macs
macs2 callpeak -t $projPath/alignment/bam/${histName}_bowtie2.mapped.bam --outdir $projPath/peakCalling/macs -n ${histName} --broad -g hs --keep-dup all -f BAMPE -q 0.05 --tempdir tmp

# ## duplicates
# ml picard/2.18.29-Java
# picardCMD="java -jar $EBROOTPICARD/picard.jar"
# mkdir -p $projPath/alignment/rmDuplicate/picard_summary

# $picardCMD SortSam I=$projPath/alignment/sam/${histName}_bowtie2.sam O=$projPath/alignment/sam/${histName}_bowtie2.sorted.sam SORT_ORDER=coordinate
# ## mark duplicates
# $picardCMD MarkDuplicates I=$projPath/alignment/sam/${histName}_bowtie2.sorted.sam O=$projPath/alignment/rmDuplicate/${histName}_bowtie2.sorted.dupMarked.sam METRICS_FILE=$projPath/alignment/rmDuplicate/picard_summary/${histName}_picard.dupMark.txt

# ## remove duplicates
# $picardCMD MarkDuplicates I=$projPath/alignment/sam/${histName}_bowtie2.sorted.sam O=$projPath/alignment/rmDuplicate/${histName}_bowtie2.sorted.rmDup.sam REMOVE_DUPLICATES=true METRICS_FILE=$projPath/alignment/rmDuplicate/picard_summary/${histName}_picard.rmDup.txt

# ## filter by the fragment length
# awk -v lb=$lowerBound -v ub=$upperBound 'BEGIN { FS="\t"; S1=lb*lb; S2=ub*ub } /^@/ { print $0; next } { if ($9*$9 >= S1 && $9*$9 <= S2) print $0}' $projPath/alignment/sam/${histName}_bowtie2.sam >$projPath/alignment/sam/${histName}_bowtie2.fragLenFilter.sam


## deeptools heatmap visualization
# ml SAMtools/1.10-foss-2016b
# ml deepTools/3.3.1-foss-2016b-Python-3.7.4

## convert into bigwig files
# mkdir -p $projPath/alignment/bigwig
# samtools sort -o $projPath/alignment/bam/${histName}.sorted.bam $projPath/alignment/bam/${histName}_bowtie2.mapped.bam
# samtools index $projPath/alignment/bam/${histName}.sorted.bam
# bamCoverage -b $projPath/alignment/bam/${histName}.sorted.bam -o $projPath/alignment/bigwig/${histName}_raw.bw

## Gene
# computeMatrix scale-regions -S $projPath/alignment/bigwig/K27me3_rep1_raw.bw \
#                                $projPath/alignment/bigwig/K27me3_rep2_raw.bw \
#                                $projPath/alignment/bigwig/K4me3_rep1_raw.bw \
#                                $projPath/alignment/bigwig/K4me3_rep2_raw.bw \
#                               -R $projPath/data/hg38_gene/hg38_gene.tsv \
#                               --beforeRegionStartLength 3000 \
#                               --regionBodyLength 5000 \
#                               --afterRegionStartLength 3000 \
#                               --skipZeros -o $projPath/data/hg38_gene/matrix_gene.mat.gz

# plotHeatmap -m $projPath/data/hg38_gene/matrix_gene.mat.gz -out $projPath/data/hg38_gene/Histone_gene.png --sortUsing sum

## Peak
# histName=$1
# repName=$2
# awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' $projPath/peakCalling/SEACR/${histName}_${repName}_seacr_control.peaks.stringent.bed >$projPath/peakCalling/SEACR/${histName}_${repName}_seacr_control.peaks.summitRegion.bed

# computeMatrix reference-point -S $projPath/alignment/bigwig/${histName}_${repName}_raw.bw \
# 	      -R $projPath/peakCalling/SEACR/${histName}_${repName}_seacr_control.peaks.summitRegion.bed \
#               --skipZeros -o $projPath/peakCalling/SEACR/${histName}_${repName}_SEACR.mat.gz -p 4 -a 3000 -b 3000 --referencePoint center

# plotHeatmap -m $projPath/peakCalling/SEACR/${histName}_${repName}_SEACR.mat.gz -out $projPath/peakCalling/SEACR/${histName}_${repName}_SEACR_heatmap.png --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${histName} ${repName}"
