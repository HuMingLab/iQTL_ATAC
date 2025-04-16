#!/bin/bash

mouse=$1
NAME=$4
CC1=$2
CC2=$3
PAT=$5
DATA=
OUT=
genome=
T=88

F1=${DATA}/${PAT}_R1_001.fastq.gz
F2=${DATA}/${PAT}_R2_001.fastq.gz

if [ ! -d "${OUT}/1.clean" ]; then
  mkdir -p "${OUT}/1.clean"
  echo "Directory created: ${OUT}/1.clean"
else
  echo "Directory already exists: ${OUT}/1.clean"
fi
if [ ! -d "${OUT}/2.statics" ]; then
  mkdir -p "${OUT}/2.statics"
  echo "Directory created: ${OUT}/2.statics"
else
  echo "Directory already exists: ${OUT}/2.statics"
fi
if [ ! -d "${OUT}/3.bam" ]; then
  mkdir -p "${OUT}/3.bam"
  echo "Directory created: ${OUT}/3.bam"
else
  echo "Directory already exists: ${OUT}/3.bam"
fi

cd ${OUT}/3.bam/
echo -e "raw reads\t$(zcat ${F1} | wc -l | awk '{print $1/4}')" > ${OUT}/2.statics/${NAME}_unique_file.txt

#conda install trim-galore(TrimGalore-0.6.10)
trim_galore -j ${T} -q 36 --length 36  --max_n 3  --stringency 3 --phred33 --paired -a 'CTGTCTCTTATACACATCT' --fastqc  -o ${OUT}/1.clean/ ${F1} ${F2} 
mv ${OUT}/1.clean/${PAT}_R1_001_val_1.fq.gz ${OUT}/1.clean/${NAME}_R1.fastq.gz
mv ${OUT}/1.clean/${PAT}_R2_001_val_2.fq.gz ${OUT}/1.clean/${NAME}_R2.fastq.gz
echo -e "after QC (clean) reads\t$(zcat ${OUT}/1.clean/${NAME}_R1.fastq.gz | wc -l | awk '{print $1/4}')" >> ${OUT}/2.statics/${NAME}_unique_file.txt


#conda install bowtie2
if ! ls ${genome}/${CC1}x${CC2}_genome.fa.*.bt2l  1> /dev/null 2>&1; then
    bowtie2-build --threads ${T} ${genome}/${CC1}x${CC2}_genome.fa ${genome}/${CC1}x${CC2}_genome.fa
else
    echo "Bowtie2 index exists."
fi

F1a=${OUT}/1.clean/${NAME}_R1.fastq.gz
F2a=${OUT}/1.clean/${NAME}_R2.fastq.gz

echo "Starting alignment..."
bowtie2 --very-sensitive -X 2000 -k 10 -p ${T} -x ${genome}/${CC1}x${CC2}_genome.fa -1 ${F1a} -2 ${F2a} | samtools view -@ ${T} -bS - > ${NAME}_${CC1}x${CC2}_genome.bam
echo "Alignment finished!!!"

samtools sort -@ ${T} -o ${NAME}_${CC1}x${CC2}_genome_sorted.bam ${NAME}_${CC1}x${CC2}_genome.bam

#conda install bioconda::picard (https://anaconda.org/bioconda/picard)
## PICARD options. You may use picard.jar with java or install picard via bioconda. comment/uncomment the appropriate method 

# java option
java -jar picard_2.27.5.jar MarkDuplicates -I ${NAME}_${CC1}x${CC2}_genome_sorted.bam \
 -O ${NAME}_${CC1}x${CC2}_genome_sorted.dedup.bam -M  ${OUT}/2.statics/${NAME}_${CC1}x${CC2}_genome.dedup_metrics -REMOVE_DUPLICATES true 

#bioconda option

#picard MarkDuplicates -I ${NAME}_${CC1}x${CC2}_genome_sorted.bam \
# -O ${NAME}_${CC1}x${CC2}_genome_sorted.dedup.bam -M  ${OUT}/2.statics/${NAME}_${CC1}x${CC2}_genome.dedup_metrics -REMOVE_DUPLICATES true 

samtools sort -@ ${T} -n -o ${NAME}_${CC1}x${CC2}_genome_srt.bam ${NAME}_${CC1}x${CC2}_genome_sorted.dedup.bam 
#conda install pysam
echo "Starting split..."
python 0.split_unique.py ${NAME}_${CC1}x${CC2}_genome_srt.bam ${NAME}_${CC1} ${NAME}_${CC2}
cat  ${NAME}_Sstats.txt| awk 'NR==1 {print $0}' >> ${OUT}/2.statics/${NAME}_unique_file.txt
cat  ${NAME}_Sstats.txt| awk 'NR==5 {print $0}' >> ${OUT}/2.statics/${NAME}_unique_file.txt
cat  ${NAME}_Sstats.txt| awk 'NR==6 {print $0}' >> ${OUT}/2.statics/${NAME}_unique_file.txt
cat  ${NAME}_Sstats.txt| awk 'NR==7 {print $0}' >> ${OUT}/2.statics/${NAME}_unique_file.txt
cat  ${NAME}_Sstats.txt| awk 'NR==8 {print $0}' >> ${OUT}/2.statics/${NAME}_unique_file.txt

samtools sort ${OUT}/3.bam/${NAME}_${CC1}_common.bam -o ${OUT}/3.bam/sorted_${NAME}_${CC1}_common.bam
samtools sort ${OUT}/3.bam/${NAME}_${CC2}_common.bam -o ${OUT}/3.bam/sorted_${NAME}_${CC2}_common.bam
samtools sort ${OUT}/3.bam/${NAME}_${CC1}_unique.bam -o ${OUT}/3.bam/sorted_${NAME}_${CC1}_unique.bam
samtools sort ${OUT}/3.bam/${NAME}_${CC2}_unique.bam -o ${OUT}/3.bam/sorted_${NAME}_${CC2}_unique.bam

samtools index ${OUT}/3.bam/sorted_${NAME}_${CC1}_common.bam
samtools index ${OUT}/3.bam/sorted_${NAME}_${CC2}_common.bam
samtools index ${OUT}/3.bam/sorted_${NAME}_${CC1}_unique.bam
samtools index ${OUT}/3.bam/sorted_${NAME}_${CC2}_unique.bam

cd WORKING DIRECTORY

./run_extract.sh ${CC1} ${CC2} ${NAME}
./run_counts.sh ${CC1} ${CC2} ${NAME}
