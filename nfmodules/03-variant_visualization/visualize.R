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

output_1ktsv <- str_replace( string = input_tsv,
                             pattern = ".tsv",
                             replacement = ".1kgFreqs.tsv" )

output_haplodots <- str_replace( string = input_tsv,
                                 pattern = ".tsv",
                                 replacement = ".haplodots.svg" )

# read the data
variants <- vroom( file = input_tsv ) %>%
  filter( AC != 1 )

# reorder data to add index
reordered <- variants %>% 
  separate( data = .,
            col = Location,
            into = c("chrom", "pos"), remove = FALSE,
            sep = ":" ) %>% 
  mutate( pos = as.numeric( pos ) ) %>% 
  arrange( chrom, pos ) %>% 
  mutate( plot_index = 1:n( ),
          .before = Uploaded_variation )

# create barplot
forbars_data <- reordered %>%
  select( plot_index, chrom, pos, Uploaded_variation, Existing_variation, AFR_AF:SAS_AF, AF ) %>%
  pivot_longer( cols = 6:ncol(.),
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

# prepare axis
xbreaks <- reordered %>% 
  pull( plot_index )

xlabs <- reordered %>% 
  pull( Existing_variation )

# plot data
bars1k <- ggplot( data = forbars_data,
                  mapping = aes( x = plot_index,
                                 y = freq,
                                 fill = pop ) ) +
  geom_col( position = "dodge",
            mapping = aes( alpha = freq ) ) +
  geom_col( data = forbars_tmp, color = "black", position = "dodge", size = 0.2 ) +
  scale_y_continuous( labels = percent ) +
  scale_x_continuous( breaks = xbreaks,
                      labels = xlabs ) +
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

# save table for 1kGP
reordered %>% 
  select( Uploaded_variation, chrom, pos,
          Allele, Existing_variation,
          AC, AN, AF, AFR_AF:SAS_AF ) %>% 
  write.table( x = ., file = output_1ktsv, append = FALSE, quote = FALSE,
               sep = "\t", row.names = FALSE, col.names = TRUE )


# intentamos un plot de haplotipo
forhaplopaint <- reordered %>% 
  select( pos, Allele, 26:ncol(.) ) %>% 
  pivot_longer( data = .,
                cols = 3:ncol(.),
                names_to = "ID",
                values_to = "gt" ) %>% 
  mutate( gt = ifelse( test = gt == ".",
                       yes = NA,
                       no = gt ),
          homhet = case_when( str_detect( string = gt, pattern = "[[:lower:]]") ~ "heterozygous",
                              str_detect( string = gt, pattern = "[[:upper:]]") ~ "homozygous") ) %>% 
  mutate( nt = toupper( gt ) )

# create gt scale
# mygtcolors <- c( "a" = "navyblue",
#                  "A" = "navyblue",
#                  "t" = "tomato",
#                  "T" = "tomato",
#                  "c" = "gold4",
#                  "C" = "gold4",
#                  "g" = "purple",
#                  "G" = "purple" )
mygtcolors <- c( 
  "A" = "navyblue",
  
  "T" = "tomato",
  
  "C" = "gold4",
  
  "G" = "purple" )

# create a gt size scale
homsize <- 3
hetzize <- 2

# # create gt scale
# mygtsize <- c( "a" = hetszize,
#                "A" = homsize,
#                "t" = hetszize,
#                "T" = homsize,
#                "c" = hetszize,
#                "C" = homsize,
#                "g" = hetszize,
#                "G" = homsize )
mygtsize <- c( "heterozygous" = hetzize,
               "homozygous" = homsize )

# create gt scale for shapes
homshape <- 16
hetshape <- 1

# mygtshape <- c( "a" = hetshape,
#                 "A" = homshape,
#                 "t" = hetshape,
#                 "T" = homshape,
#                 "c" = hetshape,
#                 "C" = homshape,
#                 "g" = hetshape,
#                 "G" = homshape )

mygtshape <- c( "heterozygous" = hetshape,
                "homozygous" = homshape )

# create function for prettynum in x axis
myprettynum <- function( the_number ) {
  
  prettyNum( x = the_number, big.mark = "," )
  
}

# create plot
haplodots <- ggplot( data = forhaplopaint) +
  geom_segment( mapping = aes( x = min(pos),
                               xend = max(pos),
                               y = ID,
                               yend = ID ), 
                size = 0.05, lty = "dotted",
                alpha = 0.1 ) +
  geom_point( data = filter( forhaplopaint, !is.na( gt ) ),
              mapping = aes( x = pos,
                             y = ID,
                             color = nt,
                             shape = homhet,
                             size = homhet 
              ) ) +
  scale_color_manual( values = mygtcolors ) +
  scale_shape_manual( values = mygtshape ) +
  scale_size_manual( values = mygtsize ) +
  scale_x_continuous( labels = myprettynum ) +
  labs( title = paste( "SNPs in region:", titulo ) ,
        subtitle = paste( cromosoma, "from", startpos, "nt to", endpos, "nt" ),
        caption = "Removed singletons from dataset",
        x = "position",
        y = "ID",
        size = "zygosity",
        shape = "zygosity" ) +
  theme_classic( ) +
  theme( plot.title = element_text( hjust = 0.5, size = 15 ),
         plot.subtitle = element_text( hjust = 0.5, size = 12 ),
         plot.caption = element_text( size = 12 ),
         panel.grid.major.x = element_blank( ),
         panel.grid.minor.x = element_blank( ),
         # axis.text.x = element_text( angle = 90, hjust = 1, vjust = 0.5 ),
         legend.position = "right",
         legend.background = element_rect( color = "black" ) )  

#
# save as svg
ggsave( filename = output_haplodots,
        plot = haplodots,
        width = 15,
        height = 10 )