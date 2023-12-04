Noticed that there was some strange annotation in Tan2018 paper with zonal anno: different loci assigned the same gene, and gene assignment not matching UCSC or any other annotation I can find. I wanted to use the Tan loci (since transcription over those loci were used to annotate zones), but I wanted the names to be compatible with all the RNAseq I did. Gerenrally I noticed that the Tan annotation were actually better than the soria in regards to the TSS, so I feel comfortable using that annotation.

I used bedtools intersect between our soria Olfr annotation (bed file with gene name) with the Tan 2018 zonal anno, with gene name and loci. I grabbed the Tan loci, wiht the Soria names. There were a few instances where there was no overlap with the Soria file, in which case I grabbed the Tan name. Also in a few instances this resulted in duplicates. I fixed these manually on a case by case basis.

Here is the command:
#cat zonal_anno_2020_fine_INTERSECT_OR-xscript-Soria+UCSC.mm10.tsv  | awk '{FS="\t"}{OFS="\t"}{if ($11 == ".") print $1,$2,$3,$4,"0",$6,$7; else print $1,$2,$3,$11,"12",$6,$7}'| awk '{FS="\t"}{OFS="\t"}{print $4}' > zonal_anno_2020_fine_TanLocus_SoriaName.bed

To determine duplicates I ran:
#cat zonal_anno_2020_fine_INTERSECT_OR-xscript-Soria+UCSC.mm10.tsv  | awk '{FS="\t"}{OFS="\t"}{if ($11 == ".") print $1,$2,$3,$4,"0",$6,$7; else print $1,$2,$3,$11,"12",$6,$7}'| awk '{FS="\t"}{OFS="\t"}{print $4}' | uniq -d

These were the duplicate Olfrs:
Olfr432
Olfr1029
Olfr78
Olfr467
Olfr399
Olfr116

not really all of them:
Olfr279
Olfr181
Olfr362
Olfr1025-ps1
Olfr1034
Olfr1195
Olfr1263
Olfr577
Olfr18

Olfr1013
Olfr1294
Olfr1396
Olfr558
Olfr615
Olfr776
Olfr777
