/* Load VCF file */
input_testdatatsv = Channel.fromPath( "./test/data/*.tsv" )
input_tsv = Channel.fromPath( "./*.tsv" )

/* gather inputs 2 */
input_tsv
  .mix( input_testdatatsv )
  // .view( )
  .set{ module_input_tsv }

/* Read mkfile module files */
Channel
	.fromPath("./*.R")
	.toList()
	.set{ rfiles }

/* extract the regions individually */
process visualize_variants {

  input:
  file tsv from module_input_tsv
  file rscripts from rfiles

  output:
  file "*.svg"

  """
  # Viz data in R
  Rscript --vanilla visualize.R $tsv
  """
}
