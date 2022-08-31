/* Load VCF file */
input_vcf = Channel.fromPath( params.input_vcf )
input_index = Channel.fromPath( params.input_index )
input_regions = Channel.fromPath( params.input_regions )

/* gather inputs 2 */
input_vcf
  .combine( input_index )
  // .view( )
  .set{ module_input_vcf }

/* extract the regions individually */
process sep_regions {

    input:
    file regions from input_regions

    output:
    file "region_*" into results_sep_regions mode flatten

    """
    split --lines=1 --suffix-length=6 -d $regions region_
    """
}

/* gather [vcf and index] + regions */
module_input_vcf
  .combine( results_sep_regions )
  // .view( )
  .set{ module_input_forfilter }

/* Filter VCF by each region */
process filter_vcf {

  input:
  file vcf_set from module_input_forfilter

  output:
  file "*.filtered.vcf" into results_filter_vcf

  """
  regionposition=\$( cut -d" " -f1 region_* )
  regionname=\$( cut -d" " -f2 region_* )
  bcftools view *.vcf.gz \$regionposition \
  | bcftools norm --multiallelics -any \
  | bcftools view --types snps \
  | bcftools annotate \
    --remove FORMAT,QUAL,^INFO/AC,INFO/AF,INFO/AN \
    > \$regionname.filtered.vcf
  """
}
