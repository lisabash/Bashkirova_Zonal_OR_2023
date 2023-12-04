#!/bin/bash

# used to make heatmaps of ChIP signal over OR genes
# same parameters are used to make plots in Figure2, Supp Figure2, Figure 6 and Supp Figure 6
# input is bigwig files merged from two replicates (generated with HOMER) and bed regions for OR gene bodies

computeMatrix scale-regions \
	      -S    mergedReps.nchip.Z1_OSN_tk79.pileup.bigWig mergedReps.nchip.Z2_OSN_tk79.pileup.bigWig mergedReps.nchip.Z5_OSN_tk79.pileup.bigWig  \
	      -R zonal_anno_classI_bed6.bed zonal_anno_zone1_bed6.bed zonal_anno_zone2_bed6.bed zonal_anno_zone3_bed6.bed zonal_anno_zone4_bed6.bed zonal_anno_zone5_bed6.bed \
	      --regionBodyLength 6000 \
	      --beforeRegionStartLength 2000 \
	      --afterRegionStartLength 2000 \
	      --missingDataAsZero \
	      --binSize 10 \
	      -out matrix.zonalOSN.tk79.merge_bed6.bin-10.flank2000.mat

plotHeatmap matrix.zonalOSN.tk79.merge_bed6.bin-10.flank2000.mat  --sortUsing mean --yAxisLabel "coverage" --samplesLabel z1OSN_merge z2OSN_merge z5OSN_merge  --legendLocation upper-right --colorMap "Reds" --heatmapHeight 12 \
	    --outFileSortedRegions deeptools.zonalOSN.tk79.merge_bed6.bin-10.flank2000.bed \
	    --outFileName heatmap.zonalOSN.tk79.merge_bed6.bin-10.flank2000_height12.pdf

