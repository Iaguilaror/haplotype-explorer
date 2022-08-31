# read libs
library( "dplyr" )
library( "vcfR" )
library( "tidyr" )
library( "stringr" )
library( "ggplot2" )
library( "vroom" )

## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## Uncomment For debugging only
## Comment for production mode only
# args[1] <- "test/data/FADS1-ENSG00000149485.filtered.vcf" ## unnanotated vcf with Genotypes AN AC AF#
# args[2] <- "test/results/FADS1-ENSG00000149485.filtered.vcf.vep.annotated.tsv" ## anotated TSV from vep#

# pass to named args
input_vcf <- args[1]
input_tsv <- args[2]

# rename the output
output_tsv <- str_replace( string = input_tsv,
                           pattern = ".tsv",
                           replacement = ".recoded.tsv" )

#read vcf
thevcf <- read.vcfR( input_vcf )

# pass info data to DF
thevariants <- getFIX( thevcf, getINFO = TRUE )[ , c( 1,2,4,5,8) ] %>% 
  as.data.frame( ) %>% 
  mutate( Uploaded_variation = paste( CHROM, POS, REF, sep = "_" ),
          Uploaded_variation = paste( Uploaded_variation, ALT, sep = "/" ) ) %>% 
  select( ALT:Uploaded_variation ) %>% 
  separate( data = ., col = INFO, into = c( "AC", "AF", "AN" ), sep = ";"  ) %>% 
  mutate( AC = str_remove( string = AC, pattern = ".*=" ),
          AN = str_remove( string = AN, pattern = ".*=" ),
          AF = str_remove( string = AF, pattern = ".*=" ) ) %>% 
  mutate( AC = as.numeric( AC ),
          AN = as.numeric( AN ),
          AF = as.numeric( AF ) )

# pass GT data to DF
thegt <- thevcf@gt[ , ]%>% 
  as.data.frame( ) %>% 
  select( -FORMAT )

# bind cols
binded <- bind_cols(  thevariants, thegt )

# recode to letter althom( ALT uppercaps ) althet ( alt lowercaps ), refhom( . )
recoded <- binded %>% 
  pivot_longer( cols = 6:ncol(.),
                names_to = "sample",
                values_to = "gt" ) %>% 
  mutate( recode = case_when( gt == "0/1" ~ tolower( ALT ),
                              gt == "1/1" ~ toupper( ALT ),
                              gt == "0/0" ~ "." ) ) %>% 
  select( -gt ) %>% 
  pivot_wider( data = ., id_cols = 1:5,
               names_from = "sample", values_from = "recode" ) %>% 
  select( -ALT )

# Delete the vcf from memory
rm( thevcf )
rm( thegt )
rm( thevariants )
rm( binded )
gc( )

# Read the annotated file
theannotation <- vroom( input_tsv, na = "-" )

# join annotation
joined <- left_join( x = theannotation,
                     y = recoded,
                     by = "Uploaded_variation" )

# clean env
rm( recoded )
rm( theannotation )
gc( )

# write output
write.table( x = joined,
             file = output_tsv,
             append = FALSE, quote = FALSE,
             sep = "\t",
             row.names = FALSE, col.names = TRUE )

#FIN
