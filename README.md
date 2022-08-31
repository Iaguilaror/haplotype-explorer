# **haplotype-explorer**
Nextllow Pipeline to visualize variants and frequencies in a set of genome coordinates, from a VCF file

------------------------------------------------------------------------

### Workflow overview

![General Workflow](dev_notes/placeholder.drawio.png)

------------------------------------------------------------------------

### Input

* A VCF file with genotype data (i. e. sample columns)

* A tab separated file with coordinates and names of the genomic regions

* Output directory path

### Output

* Main output: many images with visualizations of variants in the requested genomic regions

------------------------------------------------------------------------

### Requirements

Compatible OS*:

* Ubuntu 20.04 LTS

### Software

| Requierment | Version | Required Commands * |
|-------------|---------|---------------------|
|conda| 4.12 |conda |
|Nextflow | 22.04.5 | nextflow run |
|VEP | 107 | nextflow run |
|R | v 4.2.1 | Rscript |

* These commands must be accessible from your $PATH (i.e. you should be able to invoke them from your command line).
** Plan9 port builds many binaries, but you ONLY need the mk utility to be accessible from your command line.

#### R packages requierments

| R package | Function |
|-----------|----------|
| dplyr |magrittr, filter(), mutate(), select(), arrange() |
| vcfR |read.vcf() |
| tidyr |separate() |
| stringr |str_remove() |
| scales | colour pallettes |
| ggplot2 | geom_col( ), ggplot() |
| vroom | vroom () |


TODO  ### Installation

Download fastq_QC_processing from Github repository:

    git clone https://github.com/Iaguilaror/fastq_QC_processing.git

------------------------------------------------------------------------

TODO  #### Test

To test fastq_QC_processing execution using test data, run:

    ./runtest.sh

Your console should print the Nextflow log for the run, once every
process has been submitted, the following message will appear:

     ======
     Basic pipeline TEST SUCCESSFUL
     ======

results for test data should be in the following directory:

    ./fastq_QC_processing/test/results/fastq_quality_control-results

------------------------------------------------------------------------

TODO  ### Usage

To run with your own data go to the pipeline directory and execute:

    nextflow run fastq_quality_control.nf \
	--input_dir $your_input_dir \
	--output_dir $your_output_directory \
	-resume \
	-with-report $your_output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $your_output_directory/`date +%Y%m%d_%H%M%S`.DAG.html  


TODO ### Important

Params for trimming must be modified in the nextflow.config file

* Meaning: trim_avgqual, trim_leading, trim_trailing, etcetera.

------------------------------------------------------------------------

For information about options and parameters, run:

    nextflow run fastq_quality_control.nf --help

------------------------------------------------------------------------

#### Authors

Israel Aguilar Ordoñez
