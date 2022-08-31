/* Load VCF file */
input_testdatavcf = Channel.fromPath( "./test/data/*.vcf" )
input_vcf = Channel.fromPath( "./*.vcf" )

/* gather inputs 2 */
input_vcf
  .mix( input_testdatavcf )
  // .view( )
  .set{ module_input_vcf }

/* Read mkfile module files */
Channel
	.fromPath("./*.R")
	.toList()
	.set{ rfiles }

/* extract the regions individually */
process vep_annotate {

    input:
    file vcf from module_input_vcf
    file rscripts from rfiles

    output:
    file "*.recoded.tsv"

    """
    # main vep annotation
    vep \
		--input_file $vcf \
		--format "vcf" \
		--output_file $vcf".vep.annotated.tsv" \
		--tab \
    --fields "Uploaded_variation,Location,Allele,Existing_variation,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,gnomADg_AFR_AF,gnomADg_AMI_AF,gnomADg_AMR_AF,gnomADg_ASJ_AF,gnomADg_EAS_AF,gnomADg_FIN_AF,gnomADg_MID_AF,gnomADg_NFE_AF,gnomADg_OTH_AF,gnomADg_SAS_AF" \
		--force_overwrite \
		--stats_file $vcf".vep.stats.html" \
		--warning_file $vcf".vep.err.txt" \
		--fork 1 \
		--species "homo_sapiens" \
    --merged \
		--assembly GRCh38 \
		--cache \
		--offline \
		--buffer_size 10000 \
    --use_given_ref \
    --pick \
    --af_1kg \
    --af_gnomadg \
    && grep -v "^##" $vcf".vep.annotated.tsv" \
    | sed 's/^#//' > tmp.build \
    && mv tmp.build $vcf".vep.annotated.tsv"
    # merge data in R (this probably should be its own process)
    Rscript --vanilla merge.R $vcf $vcf".vep.annotated.tsv"
    """
}
//
// /* gather [vcf and index] + regions */
// module_input_vcf
//   .combine( results_sep_regions )
//   // .view( )
//   .set{ module_input_forfilter }
//
// /* Filter VCF by each region */
// process filter_vcf {
//
//   input:
//   file vcf_set from module_input_forfilter
//
//   output:
//   file "*.filtered.vcf" into results_filter_vcf
//
//   """
//   regionposition=\$( cut -d" " -f1 region_* )
//   regionname=\$( cut -d" " -f2 region_* )
//   bcftools view *.vcf.gz \$regionposition \
//   | bcftools norm --multiallelics -any \
//   | bcftools view --types snps \
//   | bcftools annotate \
//     --remove FORMAT,QUAL,^INFO/AC,INFO/AF,INFO/AN \
//     > \$regionname.filtered.vcf
//   """
// }
