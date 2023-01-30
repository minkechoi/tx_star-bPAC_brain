#!/bin/bash

module purge
module load fastp
mkdir fsp.d120

find ./d120 -type f| grep _R1_ |sort > d120.read1
find ./d120 -type f| grep _R2_ |sort > d120.read2

sed 's/fastq.gz/fsp.fastq.gz/g' d13.read1 |sort|awk -F"/" '{print "./fsp.d120/"$4}' > d120.out1
sed 's/fastq.gz/fsp.fastq.gz/g' d13.read2 |sort|awk -F"/" '{print "./fsp.d120/"$4}' > d120.out2

mapfile R1 < d120.read1
mapfile R2 < d120.read2
mapfile O1 < d120.out1
mapfile O2 < d120.out2


for i in {0..14}
do
	fastp -g -i ${R1[i]/[[:space:]]/} -I ${R2[i]/[[:space:]]/} -o ${O1[i]/[[:space:]]/} -O ${O2[i]/[[:space:]]/}
done

#mapping_HISAT2
module purge
module load HISAT2

mkdir ./sam/d120
mkdir ./bam/d120
mkdir ./logs
mkdir ./logs/hisat2

find ./fsp.d120/ -type f| grep _R1_ |sort > d120.fsp.read1
find ./fsp.d120/ -type f| grep _R2_ |sort > d120.fsp.read2

mapfile R1 < d120.fsp.read1
mapfile R2 < d120.fsp.read2

awk -F"/" '{print $3}' d120.read1 |sort  > d120.ident

mapfile ids < d120.ident

for i in {0..14}
do
    hisat2 -p 16 --dta --tmo -x ../../../ref_genome/danio/fasta/danio11_107 -1 ${R1[i]/[[:space:]]/} -2 ${R2[i]/[[:space:]]/} -S ./sam/d120/${ids[i]/[[:space:]]/}.sam >& ./logs/hisat2/d120.hisat.log.$i
    
done

#sam to bam
module purge
module load SAMtools
for i in {0..14}
do
samtools sort -@ 16 -o ./bam/d120/${ids[i]/[[:space:]]/}.bam ./sam/d120/${ids[i]/[[:space:]]/}.sam
    rm ./sam/d120/${ids[i]/[[:space:]]/}.sam
    
done

#remove duplicates
module purge
module load picard

mkdir ./logs/picard_md

find ./bam/d120/ -type f| grep .bam |sort > d120.bamfiles
awk -F"/" '{print $4}' d120.bamfiles | sed 's/.bam/.md.bam/gi' |sort > d120.md_bam


mapfile inbam < d120.bamfiles
mapfile outbam < d120.md_bam


for i in {0..14}
do
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${inbam[i]/[[:space:]]/} O=./bam/d120/${outbam[i]/[[:space:]]/} VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true M=./logs/picard_md/d120.marked_dup_metrics$i.txt

done



#quantification

module purge 
module load StringTie

#merge

mkdir d120.string
#stringtie --merge -p 16 -G ../ref_genome/danio/annotation/Danio_rerio.GRCz11.99.gtf -o ./tron_brain.els.merged.gtf els.filelist.txt

#module purge 
#module load GffCompare

#mkdir ./gffcom/els
#gffcompare -r ../ref_genome/danio/annotation/Danio_rerio.GRCz11.99.gtf -o ./gffcom/els/ -G ./tron_brain.els.merged.gtf

#quantification
#module purge 
#module load StringTie/1.3.3b-foss-2016b


find ./bam/d120/ -type f| grep .md.bam |sort > d120.md.bamfiles
awk -F"/" '{print $4}' d120.md.bamfiles | sed 's/.md.bam//gi' |sort > d120.md.folder


mapfile files < d120.md.bamfiles
mapfile folders < d120.md.folder

for i in {0..14}
do
	stringtie -e -B -p 16 -G ../../../ref_genome/danio/annotation/Danio_rerio.GRCz11.107.gtf -o ./d120.string/bg.${folders[i]/[[:space:]]/}/${folders[i]/[[:space:]]/}.gtf -i ${files[i]/[[:space:]]/}
done	

