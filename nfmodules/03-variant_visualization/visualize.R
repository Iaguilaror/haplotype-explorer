# read libs
library( "vroom" )
library( "dplyr" )
library( "tidyr" )
library( "ggplot2" )
library( "ggsci" )
library( "scales" )
library( "stringr" )
library( "svglite" )

# Read args from command line
args = commandArgs(trailingOnly=TRUE)

## Uncomment For debugging only
## Comment for production mode only
# args[1] <- "test/data/FADS1-ENSG00000149485.filtered.vcf.vep.annotated.recoded.tsv" ## annotated tsv with recoded genotypes #

# pass to named args
input_tsv <- args[1]

# rename the output
output_1ksvg <- str_replace( string = input_tsv,
                             pattern = ".tsv",
                             replacement = ".1kgFreqs.svg" )

# read the data
variants <- vroom( file = input_tsv ) %>%
  filter( AC != 1 )

# create barplot
forbars_data <- variants %>%
  select( Uploaded_variation, Existing_variation, AFR_AF:SAS_AF, AF ) %>%
  pivot_longer( cols = 3:ncol(.),
                names_to = "pop",
                values_to = "freq" ) %>%
  mutate( tag = ifelse( test = is.na( Existing_variation ),
                        yes = Uploaded_variation,
                        no = Existing_variation ) )

# create subdata to highlight AF
forbars_tmp <- forbars_data %>%
  mutate(  freq = ifelse( test = pop == "AF",
                          yes = freq,
                          no = NA ) )

# prepare info for labels
titulo <- ( strsplit( x = input_tsv, split = "\\." ) %>% unlist( ) )[1]
cromdata <- variants %>%
  select( Location ) %>%
  separate( col = Location, into = c("chr", "pos"), sep = ":" ) %>%
  mutate( pos = as.numeric( pos ) )

cromosoma <- unique( cromdata$chr )
startpos <- min( cromdata$pos ) %>% prettyNum( x = ., big.mark = "," )
endpos <- max( cromdata$pos ) %>% prettyNum( x = ., big.mark = "," )

# plot data
bars1k <- ggplot( data = forbars_data,
                  mapping = aes( x = tag,
                                 y = freq,
                                 fill = pop ) ) +
  geom_col( position = "dodge",
            mapping = aes( alpha = freq ) ) +
  geom_col( data = forbars_tmp, color = "black", position = "dodge", size = 0.2 ) +
  scale_y_continuous( labels = percent ) +
  scale_fill_npg( ) +
  labs( title = paste( "SNPs in region:", titulo ) ,
        subtitle = paste( cromosoma, "from", startpos, "nt to", endpos, "nt" ),
        caption = "Removed singletons from dataset",
        x = "SNP",
        y = "Allele Frequency" ) +
  guides( alpha = "none",
          fill = guide_legend( nrow = 1 ) )  +
  theme_light( ) +
  theme( plot.title = element_text( hjust = 0.5, size = 15 ),
         plot.subtitle = element_text( hjust = 0.5, size = 12 ),
         plot.caption = element_text( size = 12 ),
         panel.grid.major.x = element_blank( ),
         axis.text.x = element_text( angle = 90, hjust = 1, vjust = 0.5 ),
         legend.position = "bottom", legend.background = element_rect( color = "black" ) )

# save as svg
ggsave( filename = output_1ksvg,
        plot = bars1k,
        width = 15,
        height = 10 )
