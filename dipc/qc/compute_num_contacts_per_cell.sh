#!/bin/bash
set -u
set -e

echo "#cell	number_of_contacts	percent_duplicates	number_of_contacts_corrected"

for cell_type in z1OMP z5OMP; do
	for cell_num in $(seq 1 48); do
		fraction=$(grep "$cell_type\.$cell_num " determine_sex.txt | cut -d ' ' -f 2)
		is_female=$(echo "$fraction > 0.03" | bc -l)
		if [ $is_female == "1" ]; then
			contacts_seg_file="merged.$cell_type.$cell_num.ca.female.contacts.seg.gz"
		else
			contacts_seg_file="merged.$cell_type.$cell_num.ca.male.contacts.seg.gz"
		fi
		percent_duplicates=$(hickit -i $contacts_seg_file -o - 2>&1 1>/dev/null | grep "duplicate rate" | cut -d ' ' -f 4 | cut -d '%' -f 1)
		number_of_contacts=$(zcat $contacts_seg_file | grep -v "^#" | wc -l)
		number_of_contacts_corrected=$(echo "$number_of_contacts * ((100 - $percent_duplicates)/100)" | bc -l)
		echo "$cell_type.$cell_num	$number_of_contacts	$percent_duplicates	$number_of_contacts_corrected"
	done
done
