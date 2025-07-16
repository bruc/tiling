#!/bin/sh

# This script runs all four PacBio examples.

# The following variable results in cleaning temporary and
# intermediate files as well as compressing the results.

clean=1

datasets="2021-scAAV-CBA-eGFP 2022-ssAAV-pV-CMF-GFP 2022-ssAAV-scAAV-mix 2023-Revio-ssAAV-pAV-CMF-GFP"

for set in $datasets
do
    echo At $(date): Starting $set
    ./run_tiling_${set}.sh > run_tiling_${set}.log.1 2>&1
    if ((clean==1))
    then
	pushd ../tiling/${set} >/dev/null
	rm *.sqlite *.log
	rm -rf old
	rmdir tmp
	cd phase1
	gzip -9 *
	popd >/dev/null
    fi
done
echo At $(date): Done
