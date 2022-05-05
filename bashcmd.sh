#! /bin/bash

# --------------------------Connecting to Costa Lab HPC-----------------------------------------------------------------------------------------------------------------

# login: ssh keshav@134.130.18.27
# pwd: Aachen198921


# Transfer files local computer to HPC:
scp -r /root/alina_rnaseq/mapping_alignment/mouse_genome_m39/ keshav@134.130.18.27:/data/keshav/alina_rnaseq/mapping_alignment/
# Transfer files from HPC to local server: https

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
# also add -m 20
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
fastqc *.fastq.gz --outdir= ./qc/
for f in *.fastq.gz; do cutadapt -j 4 -a "poly A=A{20}" --quality-cutoff 20 -m 20 -u 12 -o /data/keshav/alina_rnaseq/cutadapt/trimmed_run6/run6_trimmed_$f $f 1>> /data/keshav/alina_rnaseq/reports/report_Cutadapt_run6.txt; done


## Running multiqc on these file after generating reports of fastq files

multiqc . -o ./qc/

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

# Links and command for star genome indexing

wget https://ftp.ensembl.org/pub/release-106/gtf/mus_musculus/Mus_musculus.GRCm39.106.gtf.gz  --no-check-certificate
wget https://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz --no-check-certificate

# Command for running on costa lab cluster 
STAR  --runThreadN 4 --runMode genomeGenerate --genomeDir /data/keshav/alina_rnaseq/mapping_alignment/genome_index_fromSTAR_usingm39/ --genomeFastaFiles /data/keshav/alina_rnaseq/mapping_alignment/mouse_genome_m39/Mus_musculus.GRCm39.dna.primary_assembly.fa --sjdbGTFfile /data/keshav/alina_rnaseq/mapping_alignment/mouse_genome_m39/Mus_musculus.GRCm39.106.chr.gtf --sjdbOverhang 63 --limitGenomeGenerateRAM 83476436576

# genome generation succesfully completed!!!

# Running the mapping job
# this bash script to be run from the folder that has fastq files and also needs a folder called aligned in the same folder for all the output

cd /path/to/fastq
mkdir aligned    

#!/bin/bash

# define variables
gd=/root/alina_rnaseq/mapping_alignment/genome_index_fromSTAR_usingm39/
# get our data files
# FILES=/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_758_S6_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_756_S15_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_760_S12_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_764_S7_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_768_S17_R2_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_Ctrl2_S21_R2_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_476_S20_R2_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_757_S11_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_761_S4_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_765_S1_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_768_S8_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_754_S3_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_757_S18_R2_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_762_S5_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_766_S13_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_769_S9_R1_001.fastq
#   --readFilesCommand zcat

for f in $(ls *.fastq)
do
    echo $f
    STAR --runThreadN 4 --genomeDir $gd --readFilesIn $f --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts --outFileNamePrefix ./aligned/$f.
done

echo "done!"

