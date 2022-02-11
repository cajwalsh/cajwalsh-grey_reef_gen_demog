## Script to calculate diversity metrics of populations ##

## Base directory
workingdir=/Users/cameronwalsh/Desktop/shark_nf
cd ${workingdir}

## Create directories
mkdir results/NeEstimator
mkdir results/NeEstimator/tmp

## Run R script to get .gen files for each pop and chrom./locus table for LDNe
Rscript --vanilla NeEstimator/NeEstimator_prep.R

## Move .gen files and chrom.txt to input directory for NeEstimator runs
cd results/NeEstimator
mkdir input
mv chrom.txt input
find . -name *.gen -exec mv {} ./input \;

## Remove tmp
rm -r tmp

## Set up NeEstimator input files and run for each pop
cd ${workingdir}/NeEstimator
for genepop in $(ls ${workingdir}/results/NeEstimator/input/ | grep .*gen)
do
  cat info_orig > info
  sed "s;CHANGEDIR;${workingdir};g" info > tmpfile ; mv tmpfile info
  sed "s;CHANGEINPUT;${genepop};" info > tmpfile ; mv tmpfile info
  popname=$(echo ${genepop} | sed "s;_no_monomorphs_genepop.gen;;")
  sed "s;CHANGEPOP;${popname};" info > tmpfile ; mv tmpfile info
  ./Ne2-1M i:info o:option
done
rm info
