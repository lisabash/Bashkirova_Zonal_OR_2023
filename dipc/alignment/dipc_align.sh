#!/bin/bash
set -u
set -x

##sample with the most reads: z5OMP.25 (~4 million reads)
##sample with fewer reads (use to test pipeline): z5OMP.15 (only ~1 million reads)

#run from /media/storageA/lisa/dipc/
name=$1 #example: for files dipc.z5OMP.25.R1.fastq.gz and dipc.z5OMP.25.R2.fastq.gz the name is "z5OMP.25" 
SNP_file="/media/storageA/lisa/vcf/B6_and_129P2_vs_CAST_merged.snps.chr.txt.gz"

##need to run for Nextseq/Novaseq (2 color chemistry sequencing)
echo "Running cutadapt on sample $name..."
#cutadapt --cores=8  --minimum-length 18 --nextseq-trim=20  -o dipc.$name.ca.R1.fastq.gz -p dipc.$name.ca.R2.fastq.gz  200221_sequencing_OMP/FASTQ_Generation_02_21_2020_7_50_37-210885307/dipc.$name.R1.fastq.gz 200221_sequencing_OMP/FASTQ_Generation_02_21_2020_7_50_37-210885307/dipc.$name.R2.fastq.gz > cutadapt-report.$name.txt
cutadapt --cores=8  --minimum-length 18 --nextseq-trim=20  -o dipc.$name.ca.R1.fastq.gz -p dipc.$name.ca.R2.fastq.gz  fastq_merged/dipc.$name.R1.fastq.gz fastq_merged/dipc.$name.R2.fastq.gz > cutadapt-report.$name.txt

echo "aligning sample $name"
bwa mem -5SP -t 8 /seq/mm10/mm10.fa dipc.$name.ca.R1.fastq.gz dipc.$name.ca.R2.fastq.gz | gzip > dipc.$name.ca.sam.gz
echo "finished aligning!"

rm dipc.$name.ca.R1.fastq.gz
rm dipc.$name.ca.R2.fastq.gz

echo "generating seg file"
hickit.js sam2seg -v $SNP_file dipc.$name.ca.sam.gz | hickit.js chronly -y - | gzip > $name.ca.female.contacts.seg.gz # for female
hickit.js sam2seg -v $SNP_file dipc.$name.ca.sam.gz | hickit.js chronly - | hickit.js bedflt mm10.par.bed.gz  - | gzip > $name.ca.male.contacts.seg.gz # for male
         
# resolve haplotypes via imputation
hickit -i $name.ca.female.contacts.seg.gz -o - | bgzip > $name.ca.female.contacts.pairs.gz
hickit -i $name.ca.male.contacts.seg.gz -o - | bgzip > $name.ca.male.contacts.pairs.gz

#hickit -i contacts.pairs.gz -u -o - | bgzip > impute.pairs.gz 



#determine sex of sample
echo "calculating percent chrX"
a=`zcat $name.ca.female.contacts.pairs.gz | cut -f2 | grep "chrX" | wc -l`
b=`zcat $name.ca.female.contacts.pairs.gz | wc -l`
frac=`echo "scale=4; $a / $b" | bc`
echo "$name $frac" >> determine_sex.txt 
