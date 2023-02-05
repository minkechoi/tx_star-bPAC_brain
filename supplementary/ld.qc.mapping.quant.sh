#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-T112310 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=50G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.choi@exeter.ac.uk # email address


module purge
module load fastp
mkdir fsp.LD

find ./LD -type f| grep _R1_ |sort > LD.read1
find ./LD -type f| grep _R2_ |sort > LD.read2

sed 's/fastq.gz/fsp.fastq.gz/g' LD.read1 |sort|awk -F"/" '{print "./fsp.LD/"$4}' > LD.out1
sed 's/fastq.gz/fsp.fastq.gz/g' LD.read2 |sort|awk -F"/" '{print "./fsp.LD/"$4}' > LD.out2

mapfile R1 < LD.read1
mapfile R2 < LD.read2
mapfile O1 < LD.out1
mapfile O2 < LD.out2


for i in {0..14}
do
	fastp -g -i ${R1[i]/[[:space:]]/} -I ${R2[i]/[[:space:]]/} -o ${O1[i]/[[:space:]]/} -O ${O2[i]/[[:space:]]/}
done


#mapping_HISAT2
module purge
module load HISAT2

mkdir ./sam/LD
mkdir ./bam/LD
mkdir ./logs
mkdir ./logs/hisat2

find ./fsp.LD/ -type f| grep _R1_ |sort > LD.fsp.read1
find ./fsp.LD/ -type f| grep _R2_ |sort > LD.fsp.read2

mapfile R1 < LD.fsp.read1
mapfile R2 < LD.fsp.read2

awk -F"/" '{print $3}' LD.read1 |sort  > LD.ident

mapfile ids < LD.ident

for i in {0..14}
do
    hisat2 -p 16 --dta --tmo -x ../../../ref_genome/danio/fasta/danio11_107 -1 ${R1[i]/[[:space:]]/} -2 ${R2[i]/[[:space:]]/} -S ./sam/LD/${ids[i]/[[:space:]]/}.sam >& ./logs/hisat2/LD.hisat.log.$i
    
done

#sam to bam
module purge
module load SAMtools
for i in {0..14}
do
samtools sort -@ 16 -o ./bam/LD/${ids[i]/[[:space:]]/}.bam ./sam/LD/${ids[i]/[[:space:]]/}.sam
    rm ./sam/LD/${ids[i]/[[:space:]]/}.sam
    
done

#remove duplicates
module purge
module load picard

mkdir ./logs/picard_md

find ./bam/LD/ -type f| grep .bam |sort > LD.bamfiles
awk -F"/" '{print $4}' LD.bamfiles | sed 's/.bam/.md.bam/gi' |sort > LD.md_bam


mapfile inbam < LD.bamfiles
mapfile outbam < LD.md_bam


for i in {0..14}
do
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${inbam[i]/[[:space:]]/} O=./bam/LD/${outbam[i]/[[:space:]]/} VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true M=./logs/picard_md/LD.marked_dup_metrics$i.txt

done



#quantification

module purge 
module load StringTie

#merge

mkdir LD.string
#stringtie --merge -p 16 -G ../ref_genome/danio/annotation/Danio_rerio.GRCz11.99.gtf -o ./tron_brain.els.merged.gtf els.filelist.txt

#module purge 
#module load GffCompare

#mkdir ./gffcom/els
#gffcompare -r ../ref_genome/danio/annotation/Danio_rerio.GRCz11.99.gtf -o ./gffcom/els/ -G ./tron_brain.els.merged.gtf

#quantification
#module purge 
#module load StringTie/1.3.3b-foss-2016b


find ./bam/LD/ -type f| grep .md.bam |sort > LD.md.bamfiles
awk -F"/" '{print $4}' LD.md.bamfiles | sed 's/.md.bam//gi' |sort > LD.md.folder


mapfile files < LD.md.bamfiles
mapfile folders < LD.md.folder

for i in {0..14}
do
	stringtie -e -B -p 16 -G ../../../ref_genome/danio/annotation/Danio_rerio.GRCz11.107.gtf -o ./LD.string/bg.${folders[i]/[[:space:]]/}/${folders[i]/[[:space:]]/}.gtf -i ${files[i]/[[:space:]]/}
done	

