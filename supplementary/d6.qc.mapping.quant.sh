#!/bin/bash

module purge
module load fastp
mkdir fsp.d6

find ./d6 -type f| grep _R1_ |sort > d6.read1
find ./d6 -type f| grep _R2_ |sort > d6.read2

sed 's/fastq.gz/fsp.fastq.gz/g' d6.read1 |sort|awk -F"/" '{print "./fsp.d6/"$4}' > d6.out1
sed 's/fastq.gz/fsp.fastq.gz/g' d6.read2 |sort|awk -F"/" '{print "./fsp.d6/"$4}' > d6.out2

mapfile R1 < d6.read1
mapfile R2 < d6.read2
mapfile O1 < d6.out1
mapfile O2 < d6.out2


for i in {0..14}
do
	fastp -g -i ${R1[i]/[[:space:]]/} -I ${R2[i]/[[:space:]]/} -o ${O1[i]/[[:space:]]/} -O ${O2[i]/[[:space:]]/}
done


#mapping_HISAT2
module purge
module load HISAT2

mkdir ./sam/d6
mkdir ./bam/d6
mkdir ./logs
mkdir ./logs/hisat2

find ./fsp.d6/ -type f| grep _R1_ |sort > d6.fsp.read1
find ./fsp.d6/ -type f| grep _R2_ |sort > d6.fsp.read2

mapfile R1 < d6.fsp.read1
mapfile R2 < d6.fsp.read2

awk -F"/" '{print $3}' d6.read1 |sort  > d6.ident

mapfile ids < d6.ident

for i in {0..14}
do
    hisat2 -p 16 --dta --tmo -x ../../../ref_genome/danio/fasta/danio11_107 -1 ${R1[i]/[[:space:]]/} -2 ${R2[i]/[[:space:]]/} -S ./sam/d6/${ids[i]/[[:space:]]/}.sam >& ./logs/hisat2/d6.hisat.log.$i
    
done

#sam to bam
module purge
module load SAMtools
for i in {0..14}
do
samtools sort -@ 16 -o ./bam/d6/${ids[i]/[[:space:]]/}.bam ./sam/d6/${ids[i]/[[:space:]]/}.sam
    rm ./sam/d6/${ids[i]/[[:space:]]/}.sam
    
done

#remove duplicates
module purge
module load picard

mkdir ./logs/picard_md

find ./bam/d6/ -type f| grep .bam |sort > d6.bamfiles
awk -F"/" '{print $4}' d6.bamfiles | sed 's/.bam/.md.bam/gi' |sort > d6.md_bam


mapfile inbam < d6.bamfiles
mapfile outbam < d6.md_bam


for i in {0..14}
do
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${inbam[i]/[[:space:]]/} O=./bam/d6/${outbam[i]/[[:space:]]/} VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true M=./logs/picard_md/d6.marked_dup_metrics$i.txt

done



#quantification

module purge 
module load StringTie


mkdir d6.string


find ./bam/d6/ -type f| grep .md.bam |sort > d6.md.bamfiles
awk -F"/" '{print $4}' d6.md.bamfiles | sed 's/.md.bam//gi' |sort > d6.md.folder


mapfile files < d6.md.bamfiles
mapfile folders < d6.md.folder

for i in {0..14}
do
	stringtie -e -B -p 16 -G ../../../ref_genome/danio/annotation/Danio_rerio.GRCz11.107.gtf -o ./d6.string/bg.${folders[i]/[[:space:]]/}/${folders[i]/[[:space:]]/}.gtf -i ${files[i]/[[:space:]]/}
done	

