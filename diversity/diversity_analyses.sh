## Script to calculate diversity metrics of populations ##

## Base directory
cd /Users/cameronwalsh/Desktop/shark_nf/
mkdir results/diversity

## Get projections decided upon from previous -a easySFS runs
projs=$(tail -n1 results/easySFSpreview/all_sfs_projs.csv)

## Run easySFS on entire dataset with chosen projections
easySFS.py \
-i data/Walsh_Camblyrhynchos.vcf \
-p data/vcf_pops_strata_nogb02.txt \
-o results/diversity/easySFS \
-a \
--proj=${projs}

## Change directories and rename them as required for Paolo's "Stats_from_SFS_TD.R" script
cd results/diversity/easySFS/fastsimcoal2
ls *MAFpop0.obs | sed "s;_MAFpop0.obs;;g" > pops.txt
for pop in $(cat pops.txt)
do
  tail -n2 ${pop}_MAFpop0.obs | head -n1 > ${pop}.sfs
done
rm pops.txt
mkdir ../../SFSs
mv *.sfs ../../SFSs
cd ../../SFSs

## Run Stats_from_SFS_TD.R
Rscript --vanilla ../../../diversity/Stats_from_SFS_TD.R

## Move raw diversity stats up and calculate per site values
mv Diversity_stats.txt ..
cd ..
RScript --vanilla ../../diversity/per_site_div_recalib.R

## Run both correlation analyses
RScript --vanilla ../../diversity/wide_geographic_diversity.R
RScript --vanilla ../../diversity/narrow_geographic_diversity.R
