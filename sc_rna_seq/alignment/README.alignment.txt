Realign scRNAseq experiments prepared by NYGC using our own .gtf file that contains better annotation of of the OR genes. Alignment is done using the dropseq pipeline found at: https://github.com/broadinstitute/Drop-seq/blob/master/doc/Drop-seq_Alignment_Cookbook.pdf using version 2.3.0 located at: /media/storageA/lisa/2020-scRNAseq/dropseq_v2_software/Drop-seq_tools-2.3.0/


###################
###  ALIGNMENT ####
###################

1) Generate input unaligned-queryname-sorted.bam

example:
#java -jar /seq/picard/dist/picard.jar FastqToSam \
#       F1=Plate-00103-Q1-Pool.R1.fastq.gz \
#       F2=Plate-00103-Q1-Pool.R2.fastq.gz \
#       O=unaligned_00103-Q1-Pool.bam \
#       SM=00103-Q1-Pool \

2) run dropseq pipeline.

#bash /media/storageA/lisa/2020-scRNAseq/dropseq_v2_software/Drop-seq_tools-2.3.0/Drop-seq_alignment.sh -g /media/storageA/kevin/STAR_new -r /seq/mm10/iGenome/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -d /media/storageA/lisa/2020-scRNAseq/dropseq_v2_software/Drop-seq_tools-2.3.0/ -t /media/storageA/lisa/2020-scRNAseq/realign/temp_files_pool1 -o Plate-00103-bam -s /usr/local/bin/STAR unaligned_00103-Q1-Pool.bam

3) Make metadata for tagging merged bams with exons. 

refFlat:
../dropseq_v2_software/Drop-seq_tools-2.3.0/ConvertToRefFlat ANNOTATIONS_FILE=genes_soria_pcdh.gtf SEQUENCE_DICTIONARY=/seq/mm10/iGenome/UCSC/mm10/Sequence/WholeGenomeFasta/genome.dict OUTPUT=genes_soria_pcdh.refFlat

reduced GTF:
../dropseq_v2_software/Drop-seq_tools-2.3.0/ReduceGtf GTF=genes_soria_pcdh.gtf  SEQUENCE_DICTIONARY=/seq/mm10/iGenome/UCSC/mm10/Sequence/WholeGenomeFasta/genome.dict OUTPUT=genes_soria_pcdh_reduced.gtf

intervals:
../dropseq_v2_software/Drop-seq_tools-2.3.0/CreateIntervalsFiles SEQUENCE_DICTIONARY=/seq/mm10/iGenome/UCSC/mm10/Sequence/WholeGenomeFasta/genome.dict REDUCED_GTF=genes_soria_pcdh_reduced.gtf PREFIX=genes_soria_pcdh OUTPUT=.

4) tag bam with gene annotations.

../../dropseq_v2_software/Drop-seq_tools-2.3.0/TagReadWithInterval I=merged.bam O=gene_tagged.bam TMP_DIR=. INTERVALS=../../dropseq_metadata/genes_soria_pcdh.genes.intervals TAG=XG
../../dropseq_v2_software/Drop-seq_tools-2.3.0/TagReadWithGeneFunction O=function_tagged.bam ANNOTATIONS_FILE=../../dropseq_metadata/genes_soria_pcdh.refFlat INPUT=gene_tagged.bam

5) Detect bead synthesis. 
../../dropseq_v2_software/Drop-seq_tools-2.3.0/DetectBeadSubstitutionErrors INPUT=function_tagged.bam OUTPUT=substitution_repaired.bam TMP_DIR=. MIN_UMIS_PER_CELL=20 OUTPUT_REPORT=substitution_error_report.txt

../../dropseq_v2_software/Drop-seq_tools-2.3.0/DetectBeadSynthesisErrors INPUT=substitution_repaired.bam MIN_UMIS_PER_CELL=20 OUTPUT_STATS=synthesis_error_stats.txt SUMMARY=synthesis_error_summary.txt REPORT=synthesis_error_report.txt CREATE_INDEX=true TMP_DIR=. OUTPUT=final.bam


6) q30 filtering. 
samtools view -bq 30 final.bam -o Plate103-Q1-Pool.final.q30.bam
samtools index Plate103-Q1-Pool.final.q30.bam


#######################
#FOR plates 230 and 234
#######################
Hardcoded location of refFlat and genes.invervals file so that the pipeline doesn't crash and I don't have to do anything manually.
File is named Drop-seq_alignment_Lisa.sh


#bash /media/storageA/lisa/2020-scRNAseq/dropseq_v2_software/Drop-seq_tools-2.3.0/Drop-seq_alignment_Lisa.sh -g /media/storageA/kevin/STAR_new -r /seq/mm10/iGenome/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -d /media/storageA/lisa/2020-scRNAseq/dropseq_v2_software/Drop-seq_tools-2.3.0/ -t /media/storageA/lisa/2020-scRNAseq/realign/temp_files_230_pool1 -o Plate-00230-Q1-bam -s /usr/local/bin/STAR unaligned_00230-Q1-Pool.bam
