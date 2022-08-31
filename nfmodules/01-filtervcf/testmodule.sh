# remove previous tests
rm -rf .nextflow.log* work

# remove previous results
rm -rf test/results

# create a results dir
mkdir -p test/results

# run nf script
nextflow run main.nf \
  --input_vcf test/data/FADSregion.vcf.gz \
  --input_index test/data/FADSregion.vcf.gz.tbi \
  --input_regions test/references/regions.txt

# move module results and move to test/results
mv work/*/*/*.filtered.vcf test/results/ \
&& rm -rf work                # delete workdir only if final results were found
