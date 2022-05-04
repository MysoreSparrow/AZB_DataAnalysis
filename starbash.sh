#!/bin/bash

# define variables
gd=/root/alina_rnaseq/mapping_alignment/genome_index_fromSTAR_usingm39/
# get our data files
FILES=/root/alina_rnaseq/mapping_alignment/trimmed_run5/
#   --readFilesCommand zcat

for f in $FILES
do
    echo $f
    STAR --runThreadN 4 --genomeDir $gd --readFilesIn $f --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts --outFileNamePrefix $mapped"aligned_"
done

echo "done!"
