#!/bin/bash

# HOMER script to quantify histone ChIP reads over OR gene bodies
# takes as input ChIP bigwig files and bed files with regions to compute tag count for
# used bigwigs from two merged replicates of H3K9me3 and H3K79me3 native ChIPs (Figure 2B)

bed_dir="zonal_paper/annotation"
tk79_homer="homer.mergedReps.nchip.OE_tk79"
tk9_homer="homer.mergedReps.nchip.OE_tk9"

annotatePeaks.pl $bed_dir/zonal_anno_classI.bed mm10 -size given -noann -d $tk79_homer -d $tk9_homer > classI_OE_tk79_tk9_annotation.tsv
annotatePeaks.pl $bed_dir/zonal_anno_zone1.bed mm10 -size given -noann -d $tk79_homer -d $tk9_homer > zone1_OE_tk79_tk9_annotation.tsv
annotatePeaks.pl $bed_dir/zonal_anno_zone2.bed mm10 -size given -noann -d $tk79_homer -d $tk9_homer > zone2_OE_tk79_tk9_annotation.tsv
annotatePeaks.pl $bed_dir/zonal_anno_zone3.bed mm10 -size given -noann -d $tk79_homer -d $tk9_homer > zone3_OE_tk79_tk9_annotation.tsv
annotatePeaks.pl $bed_dir/zonal_anno_zone4.bed mm10 -size given -noann -d $tk79_homer -d $tk9_homer > zone4_OE_tk79_tk9_annotation.tsv
annotatePeaks.pl $bed_dir/zonal_anno_zone5.bed mm10 -size given -noann -d $tk79_homer -d $tk9_homer > zone5_OE_tk79_tk9_annotation.tsv


cat classI_OE_tk79_tk9_annotation.tsv | sed 1d | awk '{FS="\t"}{OFS="\t"}{print $20/($4-$3), $21/($4-$3), "classI"}' > classI_OE_tk79_tk9_tagcount.tsv
cat zone1_OE_tk79_tk9_annotation.tsv | sed 1d | awk '{FS="\t"}{OFS="\t"}{print $20/($4-$3), $21/($4-$3), "zone1"}' > zone1_OE_tk79_tk9_tagcount.tsv
cat zone2_OE_tk79_tk9_annotation.tsv | sed 1d | awk '{FS="\t"}{OFS="\t"}{print $20/($4-$3), $21/($4-$3), "zone2"}' > zone2_OE_tk79_tk9_tagcount.tsv
cat zone3_OE_tk79_tk9_annotation.tsv | sed 1d | awk '{FS="\t"}{OFS="\t"}{print $20/($4-$3), $21/($4-$3), "zone3"}' > zone3_OE_tk79_tk9_tagcount.tsv
cat zone4_OE_tk79_tk9_annotation.tsv | sed 1d | awk '{FS="\t"}{OFS="\t"}{print $20/($4-$3), $21/($4-$3), "zone4"}' > zone4_OE_tk79_tk9_tagcount.tsv
cat zone5_OE_tk79_tk9_annotation.tsv | sed 1d | awk '{FS="\t"}{OFS="\t"}{print $20/($4-$3), $21/($4-$3), "zone5"}' > zone5_OE_tk79_tk9_tagcount.tsv


cat classI_OE_tk79_tk9_tagcount.tsv zone1_OE_tk79_tk9_tagcount.tsv zone2_OE_tk79_tk9_tagcount.tsv zone3_OE_tk79_tk9_tagcount.tsv zone4_OE_tk79_tk9_tagcount.tsv zone5_OE_tk79_tk9_tagcount.tsv > zones_tagdensity.txt