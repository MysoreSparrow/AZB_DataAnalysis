#!/bin/bash

# define variables
gd=/root/alina_rnaseq/mapping_alignment/genome_index_fromSTAR_usingm39/
# get our data files
FILES=/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_758_S6_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_756_S15_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_760_S12_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_764_S7_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_768_S17_R2_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_Ctrl2_S21_R2_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_476_S20_R2_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_757_S11_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_761_S4_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_765_S1_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_768_S8_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_754_S3_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_757_S18_R2_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_762_S5_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_766_S13_R1_001.fastq,/root/alina_rnaseq/mapping_alignment/trimmed_run5/run5_trimmed_769_S9_R1_001.fastq
#   --readFilesCommand zcat

for f in $FILES
do
    echo $f
    STAR --runThreadN 4 --genomeDir $gd --readFilesIn $f --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts 
done

echo "done!"
