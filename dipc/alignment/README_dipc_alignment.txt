Readme for Lisa's dipc analysis

ENVIRONMENT
Activate the conda environment containing Python 2:
conda activate py2

SOFTWARE VERSIONS
/media/storageA/lisa/dipc/hickit-0.1.1_x64-linux/hickit
/media/storageA/lisa/dipc/hickit-0.1.1_x64-linux/hickit.js


RUN ALIGNMENT
1) To align individual cells run dipc_align.sh. (Can also make a master bash script to run alignments overnight on multiple cells one after the other). This gets you through at this stage the sex of individual cells is determined and written to the "determine_sex.txt" file.

2) If the cell is female run generate_3D_female.sh, and if it is male run generate_3D_male.sh.
