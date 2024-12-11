#!/bin/bash

CC1=$1
CC2=$2
name=$3


python Extract.py --bam_file ${CC1}x${CC2}_${name}/3.bam/sorted_${name}_${CC1}_unique.bam --bed_file ../../liftover/${CC1}_liftover_ATAC_peaks.txt --chr_label  '' --out_file ${CC1}x${CC2}_${name}/${CC1}_${name}_counts.txt
python Extract.py --bam_file ${CC1}x${CC2}_${name}/3.bam/sorted_${name}_${CC2}_unique.bam --bed_file ../../liftover/${CC2}_liftover_ATAC_peaks.txt --chr_label  '' --out_file ${CC1}x${CC2}_${name}/${CC2}_${name}_counts.txt
python Extract.py --bam_file ${CC1}x${CC2}_${name}/3.bam/sorted_${name}_${CC1}_common.bam --bed_file ../../liftover/${CC1}_liftover_ATAC_peaks.txt --chr_label  '' --out_file ${CC1}x${CC2}_${name}/${CC1}_${name}_common_counts.txt
python Extract.py --bam_file ${CC1}x${CC2}_${name}/3.bam/sorted_${name}_${CC2}_common.bam --bed_file ../../liftover/${CC2}_liftover_ATAC_peaks.txt --chr_label  '' --out_file ${CC1}x${CC2}_${name}/${CC2}_${name}_common_counts.txt
