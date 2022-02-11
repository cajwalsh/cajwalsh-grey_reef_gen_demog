absolute_path = "/Users/cameronwalsh/Desktop/shark_nf"

sample_names = file("$absolute_path/data/vcf_pops_strata_nogb02.txt").readLines()

Channel
  .from(sample_names)
  .map { n -> (n - ~/.*\s/) }
  .unique()
  .set { populations }

Channel
  .fromPath( [ "$absolute_path/data/vcf_pops_strata_nogb02.txt",
  "$absolute_path/data/Walsh_Camblyrhynchos.vcf" ] )
  .collect()
  .set { data_files }

populations
  .combine(data_files)
  .set { per_pops }

process easySFSpreview {
  publishDir "$absolute_path/results/easySFSpreview/$pop", mode: 'copy'
  tag "${pop}"

  input:
  set val( pop ), file( samples ), file( vcf ) from per_pops

  output:
  set file( "IDs.txt" ), file( "SFS.vcf" ),
      file ( "SFS_preview_thin.txt" ), file ( "SFS_preview_all.txt" )

  script:
  """
  cat $samples | \
  grep $pop > IDs.txt
  cut -f1 IDs.txt > inds_to_keep.txt

  vcftools \
    --vcf $vcf \
    --keep inds_to_keep.txt \
    --recode \
    --stdout | \
  gzip > first_step.vcf.gz

  vcftools \
    --gzvcf first_step.vcf.gz \
    --mac 1 \
    --recode \
    --stdout 2> sites.log 1> SFS.vcf

  rm first_step.vcf.gz
  rm inds_to_keep.txt

  easySFS.py -i SFS.vcf -p IDs.txt --preview > SFS_preview_thin.txt
  easySFS.py -i SFS.vcf -p IDs.txt -a --preview > SFS_preview_all.txt
  """
}
