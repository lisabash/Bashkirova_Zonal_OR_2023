#!/bin/bash
set -u
set -e

echo "#cell	rms-rmsd	median-rmsd"
for dir in $(ls | grep "^3D")
do
	log_file=$(find $dir -name "*20k.align.log")
	if [ $log_file ]; then
		rms_rmsd=$(grep "RMS RMSD" $log_file  | cut -d ' ' -f 4)
		median_rmsd=$(grep "median RMSD" $log_file | cut -d ' ' -f 4)
		echo "$dir	$rms_rmsd	$median_rmsd"
	else
		echo "#WARNING did NOT find log file for cell dir $dir"
	fi
done
