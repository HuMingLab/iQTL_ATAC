#!/bin/bash

CC1=$1
CC2=$2
mouse=$3

python allelic_counts.py ${CC1}x${CC2}_${mouse}/${CC1}_${mouse}_counts.txt ${CC1}x${CC2}_${mouse}/${CC2}_${mouse}_counts.txt ${CC1}x${CC2}_${mouse}/${CC1}_${mouse}_common_counts.txt ${CC1}x${CC2}_${mouse}/${CC2}_${mouse}_common_counts.txt $CC1 $CC2 /home/mishras10/ENCODE-ATAC-60mice/clustering_analysis/unionpeaks.all.bed results/${CC1}x${CC2}_${mouse}_allelic_total_counts.txt
