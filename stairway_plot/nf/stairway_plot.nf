absolute_path = "/Users/cameronwalsh/Desktop/shark_nf"

Channel
  .fromPath("$absolute_path/results/easySFSpreview/sp_sfs_projs.csv")
  .splitCsv()
  .set { pop_projs }

Channel
  .fromPath( "$absolute_path/stairway_plot/blueprint.blueprint" )
  .set { blueprint }

pop_projs
  .combine(blueprint)
  .set { per_pops }

process blueprint_prep {
  publishDir "$absolute_path/results/stairway_plot/output/$pop", mode: 'copy'
  tag "${pop}"

  input:
  set val( pop ), val( proj ), file ( blueprint ) from per_pops

  output:
  set val( pop ), file( "final.blueprint" ) into stairway_channel
  set file ("thinned_SFS.txt"), file ("recalibration_values.txt")

  script:
  """
  easySFS.py -i ${absolute_path}/results/easySFSpreview/${pop}/SFS.vcf \
  -p ${absolute_path}/results/easySFSpreview/${pop}/IDs.txt --proj=${proj}

  cp output/fastsimcoal2/${pop}_MAFpop0.obs thinned_SFS.txt
  rm -r output

  cp $blueprint new1.blueprint

  sed "s/CHANGEPOP/$pop/g" new1.blueprint > new2.blueprint

  sed "s/CHANGESEQS/${proj}/" new2.blueprint > new3.blueprint

  Rscript --vanilla $absolute_path/stairway_plot/recalibrations_stairway.R $pop $absolute_path | \
  cut -f2-4 -d' ' | sed 's/"//g' > recalibration_values.txt
  L_pop_thinned=\$(tail -n1 recalibration_values.txt | cut -f2 | cut -f3 -d' ')
  sed "s/CHANGELENGTH/\${L_pop_thinned}/" new3.blueprint > new4.blueprint

  nsamp=\$(echo "${proj}/2" | bc)
  nsamp_plus_1=\$(echo "1 + \${nsamp}" | bc)
  SFS_text=\$(tail -n2 thinned_SFS.txt | head -n1 | cut -d ' ' -f2-\${nsamp_plus_1})
  sed "s/CHANGESFS/\${SFS_text}/" new4.blueprint > new5.blueprint

  rand4=\$(echo "${proj}-2" | bc)
  rand1=\$(echo "\${rand4}/4" | bc)
  rand2=\$(echo "\${rand4}/2" | bc)
  rand3=\$(echo "3*\${rand4}/4" | bc)
  nrand=\$(echo "\${rand1}")
  sed "s/CHANGEBREAKPOINTS/\${nrand}/" new5.blueprint > new6.blueprint

  sed "s:CHANGEPATH:$absolute_path/stairway_plot/stairway_plot_v2/stairway_plot_es:" new6.blueprint > final.blueprint

  rm new*.blueprint
  """
}

//NOTE: although calculated above, I didn't bother using rand2,3,4 because
//      they never seemed useful in all previous times I used this program;
//      if you have fewer samples, it's possible they could be useful

process stairway_plot {
  publishDir "$absolute_path/results/stairway_plot/output/$pop", mode: 'copy'
  tag "${pop}"

  input:
  set val( pop ), file( blueprint ) from stairway_channel

  output:
  set file( "stairways/${pop}.final.summary" ),
      file( "stairways/${pop}.final.summary.pdf" )
  //  file( "blueprint.plot.sh" ) if using more than one # of
  //  rand breakpoints to figure out which was final after LT

  script:
  """
  java -cp $absolute_path/stairway_plot/stairway_plot_v2/stairway_plot_es Stairbuilder $blueprint
  bash ${blueprint}.sh
  """
}
