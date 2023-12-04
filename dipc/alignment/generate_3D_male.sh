#!/bin/bash
set -u
set -x

#z5OMP.8.ca.female.chrVCF.impute.pairs.gz
#run from /media/storageA/lisa/dipc/
name=$1 #example: for files dipc.z5OMP.25.R1.fastq.gz and dipc.z5OMP.25.R2.fastq.gz the name is "z5OMP.25"

# resolve haplotypes via imputation
hickit -i $name.ca.male.contacts.pairs.gz -u -o - | bgzip > $name.ca.male.impute.pairs.gz

mkdir 3D.$name.ca.male
cd 3D.$name.ca.male

# generate 3D structures (with 5 replicates)
for rep in `seq 1 5`
do
  hickit -s${rep} -M -i ../$name.ca.male.impute.pairs.gz -Sr1m -c1 -r10m -c2 -b4m -b1m -O $name.1m.${rep}.3dg -b200k -O $name.200k.${rep}.3dg -D5 -b50k -O $name.50k.${rep}.3dg -D5 -b20k -O $name.20k.${rep}.3dg
done

#######################
# convert from hickit to dip-c formats, and remove repetitive regions from 3D structures
../dip-c-master/scripts/hickit_pairs_to_con.sh ../$name.ca.male.contacts.pairs.gz
../dip-c-master/scripts/hickit_impute_pairs_to_con.sh ../$name.ca.male.impute.pairs.gz

for rep in `seq 1 5`
do
  ../dip-c-master/scripts/hickit_3dg_to_3dg_rescale_unit.sh $name.20k.${rep}.3dg
  ../dip-c-master/dip-c clean3 -c ../$name.ca.male.impute.con.gz $name.20k.${rep}.dip-c.3dg > $name.20k.${rep}.clean.3dg # remove repetitive (contact-less) regions
done

# align replicate structures and calculate RMSD (overall value in .log file)
../dip-c-master/dip-c align -o $name.aligned.20k. $name.20k.[1-5].clean.3dg 2> $name.20k.align.log > $name.20k.align.color


##### FILES FOR PYMOL MODELING
# color by chromosome number and visualize as mmCIF (viewable with pymol)
../dip-c-master/dip-c color -n ../dip-c-master/color/mm10.chr.txt $name.20k.1.clean.3dg | ../dip-c-master/dip-c vis -c /dev/stdin $name.20k.1.clean.3dg > $name.20k.1.clean.n.cif

# color OR genes vs nonOR (binary)
../dip-c-master/dip-c color -n ../dip-c-master/color/mm10.olfactory_receptors.binary.20k.txt $name.20k.1.clean.3dg | ../dip-c-master/dip-c vis -c /dev/stdin $name.20k.1.clean.3dg > $name.20k.1.clean.OR.cif

# color OR genes vs nonOR (binary) the -M feature removes particles that are missing from the data. not sure if it is better
../dip-c-master/dip-c color -n ../dip-c-master/color/mm10.olfactory_receptors.binary.20k.txt $name.20k.1.clean.3dg | ../dip-c-master/dip-c vis -M -c /dev/stdin $name.20k.1.clean.3dg > $name.20k.1.clean.OR_M.cif

# color by cpg fequency, the -M feature removes particles that are missing from the data. may need to do for ORs
../dip-c-master/dip-c color -n ../dip-c-master/color/mm10.cpg.20k.txt $name.20k.1.clean.3dg | ../dip-c-master/dip-c vis -M -c /dev/stdin $name.20k.1.clean.3dg > $name.20k.1.clean.cgp.cif
