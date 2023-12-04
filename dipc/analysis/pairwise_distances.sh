#!/bin/bash
set -u 
set -e

# INPUTS:
# - leg files and associated name files (produced with make_leg_files.py)
# - list of good cells with their sexes (../good_cells.txt)
# - 3DG data (from ../3D.../..cleaned.1.3dg) (e.g. ../3D.merged.z1OMP.10.ca.male/merged.z1OMP.10.20k.1.clean.3dg)

# this analysis produces counts for every row (e.g. OR or enhancer)

for DISTANCE_THRESHOLD in 2.5 5 10 20; do

	echo "DISTANCE THRESHOLD: ${DISTANCE_THRESHOLD}"

while read cell_name sex; do
	echo "processing cell: $cell_name ($sex)..."

	input_3dg_file="../3D.merged.$cell_name.ca.$sex/merged.$cell_name.20k.1.clean.3dg"
	ls $input_3dg_file

	# all ORs vs all ORS
	distances_file=pairwise_distances/$cell_name.$sex.ors_vs_ors.pairwise_distances.txt
	counts_file=counts/$cell_name.$sex.orsNoGI_vs_orsNoGI.dist${DISTANCE_THRESHOLD}.counts.txt
	leg1_file=zonal_anno_2020_fine_TanLocus_SoriaName_sorted.noGIOverlap_50kb.bed.leg
	leg2_file=zonal_anno_2020_fine_TanLocus_SoriaName_sorted.noGIOverlap_50kb.bed.leg
	name1_file=zonal_anno_2020_fine_TanLocus_SoriaName_sorted.noGIOverlap_50kb.bed.name
	name2_file=zonal_anno_2020_fine_TanLocus_SoriaName_sorted.noGIOverlap_50kb.bed.name
	python ../dip-c-master/scripts/network_around.py $DISTANCE_THRESHOLD $distances_file $leg1_file $leg2_file $name1_file $name2_file > $counts_file

done < ../good_cells.txt

done
