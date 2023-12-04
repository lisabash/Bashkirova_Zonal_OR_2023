#!/bin/bash
set -u 
set -e

# script to determine density of contacts at OR alleles relative to the surounding regions
# for each cell calculates ratio of 1) mean interchromosomal contact density within 200 kb of all OR pairs and 2) the mean density in surrounding regions (within 10 Mb)


# all ORs vs all ORs
python make_asymmetric_con_file_unknown_alleles.py OR_genes.noGIOverlap_50kb.bed OR_genes.noGIOverlap_50kb.bed > ors_ors_unknown_allele_no_overlap_on_both.con


python compute_contact_density_parallel.py ors_ors_unknown_allele_no_overlap_on_both.con good_cells_1.txt contacts > .tmp2.01 &
python compute_contact_density_parallel.py ors_ors_unknown_allele_no_overlap_on_both.con good_cells_2.txt contacts > .tmp2.02 &
python compute_contact_density_parallel.py ors_ors_unknown_allele_no_overlap_on_both.con good_cells_3.txt contacts > .tmp2.03 &
python compute_contact_density_parallel.py ors_ors_unknown_allele_no_overlap_on_both.con good_cells_4.txt contacts > .tmp2.04 &
python compute_contact_density_parallel.py ors_ors_unknown_allele_no_overlap_on_both.con good_cells_5.txt contacts > .tmp2.05 &
python compute_contact_density_parallel.py ors_ors_unknown_allele_no_overlap_on_both.con good_cells_6.txt contacts > .tmp2.06 &
python compute_contact_density_parallel.py ors_ors_unknown_allele_no_overlap_on_both.con good_cells_7.txt contacts > .tmp2.07 &
python compute_contact_density_parallel.py ors_ors_unknown_allele_no_overlap_on_both.con good_cells_8.txt contacts > .tmp2.08 &
python compute_contact_density_parallel.py ors_ors_unknown_allele_no_overlap_on_both.con good_cells_9.txt contacts > .tmp2.09 &
python compute_contact_density_parallel.py ors_ors_unknown_allele_no_overlap_on_both.con good_cells_10.txt contacts > .tmp2.10 &
python compute_contact_density_parallel.py ors_ors_unknown_allele_no_overlap_on_both.con good_cells_11.txt contacts > .tmp2.11 &
python compute_contact_density_parallel.py ors_ors_unknown_allele_no_overlap_on_both.con good_cells_12.txt contacts > .tmp2.12 &
wait
echo "waiting.."
cat .tmp2.* > ors_ors_unknown_allele_contacts_densities_no_overlap_on_both.txt

