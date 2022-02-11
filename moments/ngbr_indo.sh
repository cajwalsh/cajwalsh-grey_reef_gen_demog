## Script to get 2D-SFSs and recalibrated L for ngbr_indo moments analysis ##

## Base directory
cd /Users/cameronwalsh/Desktop/shark_nf/
mkdir results/moments

## Run script to get vcfs to make 2D-SFSs
Rscript --vanilla moments/vcfs_for_joint_SFSs.R
rm results/moments/*.rad

## Make inds file
grep 'north_gbr\|indo' data/vcf_pops_strata_nogb02.txt | \
sort -r > results/moments/ngbr_indo_IDs.txt

## Run easySFS on un-thinned dataset (with projections input manually from
## easySFSpreview) to be able to calculate sumSFS, Spop, and Lpop from
## equations in Appendix S1
easySFS.py \
-i results/moments/ngbr_indo_un-thinned.vcf \
-p results/moments/ngbr_indo_IDs.txt \
-o results/moments/un-thinned_easySFS \
-a \
--proj=26,40

## Repeat easySFS on thinned dataset (with projections input manually from
## easySFSpreview) to get final thinned SFS for analysis and calculate Sthin
## and final Lpop(thinned) to recalibrate length and model timing
easySFS.py \
-i results/moments/ngbr_indo_thinned.vcf \
-p results/moments/ngbr_indo_IDs.txt \
-o results/moments/thinned_easySFS \
--proj=26,40

## This thinned SFS is the final one we will use for moments, so copy it
## to the main results directory to make it easy to find
cp results/moments/thinned_easySFS/fastsimcoal2/ngbr_indo_thinned_jointMAFpop1_0.obs \
results/moments/ngbr_indo_thinned_joint.sfs

## Now do calculations in Appendix S1 to get Lpop(thinned)
Rscript --vanilla moments/recalibrations_ngbr_indo.R | cut -f2-4 -d' ' | \
sed 's/"//g' > results/moments/recalibration_values.txt
echo "For recalibration values, check results/moments/recalibration_values.txt"
