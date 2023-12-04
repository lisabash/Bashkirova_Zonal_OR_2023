#!/bin/bash
set -e
set -u
set -x

echo "#cell	num_read_pairs	num_raw_contacts	ratio"

for cell_type in z1OMP z5OMP; do
	for cell_num in $(seq 1 48); do
		#echo "testing $cell_type.$cell_num"
		R1_fastq_file="fastq_merged/dipc.merged.$cell_type.$cell_num.R1.fastq.gz"
		contacts_seg_file="merged.$cell_type.$cell_num.ca.female.contacts.seg.gz" # using female because it doesn't matter
		if test -f $R1_fastq_file && test -f $contacts_seg_file; then
			num_raw_contacts=$(hickit --dup-dist=0 -i $contacts_seg_file -o - | grep -v "^#" | awk '{sum++}END{print sum}')
			num_read_pairs=$(zcat $R1_fastq_file | awk '{sum++}END{print sum/4}')
			ratio=$(echo "$num_raw_contacts / $num_read_pairs" | bc -l)
			echo  "$cell_type.$cell_num	$num_read_pairs	$num_raw_contacts	$ratio"
		else
			echo "#WARNING cell $cell_type.$cell_num is missing either the fastq file or contacts seg file"
		fi
	done
done
