#!/bin/bash

module purge
module load fastp
mkdir fsp.d13

find ./d13 -type f| grep _R1_ |sort > d13.read1
find ./d13 -type f| grep _R2_ |sort > d13.read2

sed 's/fastq.gz/fsp.fastq.gz/g' d13.read1 |sort|awk -F"/" '{print "./fsp.d13/"$4}' > d13.out1
sed 's/fastq.gz/fsp.fastq.gz/g' d13.read2 |sort|awk -F"/" '{print "./fsp.d13/"$4}' > d13.out2

mapfile R1 < d13.read1
mapfile R2 < d13.read2
mapfile O1 < d13.out1
mapfile O2 < d13.out2


for i in {0..14}
do
	fastp -g -i ${R1[i]/[[:space:]]/} -I ${R2[i]/[[:space:]]/} -o ${O1[i]/[[:space:]]/} -O ${O2[i]/[[:space:]]/}
done


#mapping_HISAT2
module purge
module load HISAT2

mkdir ./sam/d13
mkdir ./bam/d13
mkdir ./logs
mkdir ./logs/hisat2

find ./fsp.d13/ -type f| grep _R1_ |sort > d13.fsp.read1
find ./fsp.d13/ -type f| grep _R2_ |sort > d13.fsp.read2

mapfile R1 < d13.fsp.read1
mapfile R2 < d13.fsp.read2

awk -F"/" '{print $3}' d13.read1 |sort  > d13.ident

mapfile ids < d13.ident

for i in {0..14}
do
    hisat2 -p 16 --dta --tmo -x ../../../ref_genome/danio/fasta/danio11_107 -1 ${R1[i]/[[:space:]]/} -2 ${R2[i]/[[:space:]]/} -S ./sam/d13/${ids[i]/[[:space:]]/}.sam >& ./logs/hisat2/d13.hisat.log.$i
    
done

#sam to bam
module purge
module load SAMtools
for i in {0..14}
do
samtools sort -@ 16 -o ./bam/d13/${ids[i]/[[:space:]]/}.bam ./sam/d13/${ids[i]/[[:space:]]/}.sam
    rm ./sam/d13/${ids[i]/[[:space:]]/}.sam
    
done

#remove duplicates
module purge
module load picard

mkdir ./logs/picard_md

find ./bam/d13/ -type f| grep .bam |sort > d13.bamfiles
awk -F"/" '{print $4}' d13.bamfiles | sed 's/.bam/.md.bam/gi' |sort > d13.md_bam


mapfile inbam < d13.bamfiles
mapfile outbam < d13.md_bam


for i in {0..14}
do
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${inbam[i]/[[:space:]]/} O=./bam/d13/${outbam[i]/[[:space:]]/} VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true M=./logs/picard_md/d13.marked_dup_metrics$i.txt

done



#quantification

module purge 
module load StringTie

#merge

mkdir d13.string
#stringtie --merge -p 16 -G ../ref_genome/danio/annotation/Danio_rerio.GRCz11.99.gtf -o ./tron_brain.els.merged.gtf els.filelist.txt

#module purge 
#module load GffCompare

#mkdir ./gffcom/els
#gffcompare -r ../ref_genome/danio/annotation/Danio_rerio.GRCz11.99.gtf -o ./gffcom/els/ -G ./tron_brain.els.merged.gtf

#quantification
#module purge 
#module load StringTie/1.3.3b-foss-2016b


find ./bam/d13/ -type f| grep .md.bam |sort > d13.md.bamfiles
awk -F"/" '{print $4}' d13.md.bamfiles | sed 's/.md.bam//gi' |sort > d13.md.folder


mapfile files < d13.md.bamfiles
mapfile folders < d13.md.folder

for i in {0..14}
do
	stringtie -e -B -p 16 -G ../../../ref_genome/danio/annotation/Danio_rerio.GRCz11.107.gtf -o ./d13.string/bg.${folders[i]/[[:space:]]/}/${folders[i]/[[:space:]]/}.gtf -i ${files[i]/[[:space:]]/}
done	

