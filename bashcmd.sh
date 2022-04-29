#! /bin/bash

# --------------------------Connecting to Costa Lab HPC-----------------------------------------------------------------------------------------------------------------

# login: ssh keshav@134.130.18.27
# pwd: Aachen198921
# --------------------------Run fastqc for first time on raw sequenced files to get an idea about quality of files.-----------------------------------------------------------------------------------------------------------------

## from directory which has all s/w installed: /root/ , run the following command to run fastqc algorithm, along with specifying where the output will be stored:
## location of fastq files on wsl desktop: /root/alina_rnaseq/fastq/

## Create Output folder
mkdir /root/alina_rnaseq/fastq/output/

## Run fastqc 
#for f in *.fastq.gz; do fastqc /root/alina_rnaseq/fastq/*.fastq.gz --outdir=/root/alina_rnaseq/fastq/output/ ;done   # soemtimes this runs into a infinite lopp and thats problematic. then use the command below.
fastqc *.fastq.gz --outdir=/root/alina_rnaseq/fastq/output/

## Run multiqc
multiqc . -o /root/alina_rnaseq/fastq/output/


# --------------------------TRIMMING WITH CUTADAPT-----------------------------------------------------------------------------------------------------------------

# cutadapt to be run from /root/alina_rnaseq/fastq/

cd /root/alina_rnaseq/fastq/

for f in *.fastq.gz; do cutadapt -j 4 -a "poly A=A{20}" --quality-cutoff 20 -u 12 -o /root/alina_rnaseq/cutadapt/trimmed_run4/run4_trimmed_$f $f 1>> /root/alina_rnaseq/reports/report_Cutadapt_run4.txt; done

# Cutadapt parameters:
# - a: The sequence of the adapter is given with the -a option. You need to replace AACCGGTT with the correct adapter sequence. or also look poly A and poly T sequences. 
# - q: The -q (or --quality-cutoff) parameter can be used to trim low-quality ends from reads. If you specify a single cutoff value, the 3’ end of each read is trimmed. Used for quality trimming to remove short reads (<20nt).
# - j: To specify number of threads.
# - o: specify output directory where trimmed files are to be written by cutadapt. 

# Perform fastqc and Multiqc again on cutadapt trimmed files.

## Create Output folder
mkdir /root/alina_rnaseq/cutadapt/trimmed_run4/output/

fastqc *.fastq.gz --outdir=/root/alina_rnaseq/cutadapt/trimmed_run4/output/
# for f in *.fastq.gz; do fastqc /root/alina_rnaseq/cutadapt/trimmed_cutadapt_run2/*.fastq.gz --outdir=/root/alina_rnaseq/cutadapt/trimmed_cutadapt_run2/output/;done
# for f in *.fastq.gz; do fastqc *.fastq.gz --outdir=/root/alina_rnaseq/cutadapt/trimmed_cutadapt_run2/output/;done

# In the cluster:
fastqc *.fastq.gz --outdir=/data/keshav/alina_rnaseq/fastq/output/
for f in *.fastq.gz; do cutadapt -j 4 -a "poly A=A{20}" --quality-cutoff 20 -u 12 -o /data/keshav/alina_rnaseq/cutadapt/trimmed_run5/run5_trimmed_$f $f 1>> /data/keshav/alina_rnaseq/reports/report_Cutadapt_run5.txt; done


## Running multiqc on these file after generating reports of fastq files

multiqc . -o /root/alina_rnaseq/reports/postcutadapttrimming_run4/


# --------------------------Mapping with STAR aligner-----------------------------------------------------------------------------------------------------------------
# Basic STAR workflow consists of 2 steps:
# 1. Generating genome indexes files .
# 2. Mapping reads to the genome. 

# 1. Generating genome indexes files

# Parameters: 
# --runThreadN 4
# --runMode genomeGenerate
# --genomeDir /path/to/genomeDir
# --genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2 ...
# --sjdbGTFfile /path/to/annotations.gtf
# --sjdbOverhang 63

# nohup is used so that the process keeps going on once even if i am disconnected from the server
# & to put the process in the background


nohup STAR  --runThreadN 4 --runMode genomeGenerate --genomeDir /data/keshav/alina_rnaseq/mapping_alignment/genome_index_fromSTAR_usingm38/ --genomeFastaFiles /data/keshav/alina_rnaseq/mapping_alignment/mouse_genome_m38/Mus_musculus.GRCm38.cdna.all.fa.gz --sjdbGTFfile /data/keshav/alina_rnaseq/mapping_alignment/mouse_genome_m38/Mus_musculus.GRCm38.102.gtf.gz --sjdbOverhang 100 &



