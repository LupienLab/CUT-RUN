# CUT-RUN

1. Download peak calling tool: https://github.com/FredHutch/SEACR
2. Set tool path, project path, fastq path as following:

cores,
projPath,
fastqdir,
seacr,
trimmomaticpath,
ref,
spikeInRef,
chromSize

3. Submit the pipeline as following:

sbatch -p himem -J cutrun --export=ALL -N 1 -n 5 --mem 80G -t 1-0 --wrap "bash CUT-RUN.pipeline.sh K27me3"


Inspired from:
Zheng Y et al (2020). Protocol.io

https://yezhengstat.github.io/CUTTag_tutorial/

