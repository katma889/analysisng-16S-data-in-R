# analysisng-16S-data-in-R
analysing and visualising data in R
##### Installing required packages for microbiome data analysis with phyloseq

## phyloseq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

## Tibble
install.packages('tibble')

## GGplot2
BiocManager::install("ggplot2")

## Ape
install.packages('ape')

## Vegan
install.packages('vegan')


# Loading installed packages
library('phyloseq')
library('tibble')
library('ggplot2')
library('ape')
library('vegan')

# Setting the working directory for R analysis
setwd('/Users/mandiraka/Desktop/16S-analysis-jan2022/R-studio-all')


# Importing 
import_table <- read.table('working-dir/otu_frequency_table.tsv',header=TRUE,sep='\t',row.names=1, comment.char = "")
head(import_table)

#Converting to a matrix for Plyloseq
otumat <- as.matrix(import_table)
head(otumat)

#Creating an OTU object using the function otu_table
OTU <- otu_table(otumat, taxa_are_rows = TRUE)
head(OTU)

##Importing the taxonomy table
import_taxa <- read.table('working-dir/sintax_taxonomy.tsv',header=TRUE,sep='\t',row.names=1)
head(import_taxa)
#convert taxa table into a matrix
taxonomy <- as.matrix(import_taxa)
head(taxonomy)
#Create a taxonomy class object
TAX <- tax_table(taxonomy)
##Importing the sample metadata (where controls were removed pror import)
metadata <- read.table('working-dir/sample_metadata_weevils_only_no_control.tsv',header = T,sep='\t',row.names = 1)
metadata
tail(metadata)
# Creating a Phyloseq sample_data_class
META <- sample_data(metadata)

## Now time to create a Phyloseq object after having all three components
physeq <- phyloseq(OTU,TAX,META)
physeq
## we didnot have to convert metadata variables to ordinal which is otherwise help us for graph these fields
## Initial data inspection to check alpha rarefaction of species richness to show that sufficient sequencing has been done to detect most species
rarecurve(t(otu_table(physeq)), step=50, cex=1)
## Creating a bar plot of Abundance
plot_bar(physeq)
print(min(sample_sums(physeq)))
print(max(sample_sums(physeq)))
## Rarefied the data- 
#From the Intial look we came to know that sample A16 have at least double many reads as the other samples. So we are using rarefaction to stimulate an even number of reads per sample. Here, we will try to create rarefied version of the phyloseq object
physeq.rarefied <- rarefy_even_depth(physeq, rngseed=1, sample.size=0.9*min(sample_sums(physeq)), replace=F)
plot_bar(physeq.rarefied)
saveRDS(physeq, 'weevil_phyloseq.rds')
saveRDS(physeq.rarefied, 'weevil_phyloseq_rarefied.rds')
##saving pdf version of graph
pdf('species_richness_plot.pdf')
##setting working directory
setwd('/Users/mandiraka/Desktop/16S-analysis-jan2022/R-studio-all')
##loading the packages
library('phyloseq')
library('tibble')
library('ggplot2')
library('ape')
library('vegan')
physeq <- readRDS('weevil_phyloseq.rds')
print(physeq)
##Ploting taxonomy using plot_bar fuction so that we can plot the relative abundance of taxa across our samples
plot_bar(physeq, fill='Genus')
##plotting rarefied dataset to compare abundance across samples
plot_bar(physeq.rarefied, fill='Genus')
plot_bar(physeq.rarefied, x='Weevils', fill='Genus')
##No difference was obtained in the plot structure with the above two commands for lotting taxonomy table
plot_bar(physeq.rarefied, x='Weevils', fill='Family')
rlang::last_error()
barplot1 <- plot_bar(physeq.rarefied, fill='Order') 
barplot1 + facet_wrap(~Weevils, scales="free_x",nrow=1)
barplot1 + facet_wrap(~Weevils, scales="free_x",nrow=1)
plot_bar(physeq, fill='Phylum')
plot_bar(physeq.rarefied, fill='Phylum')
getwd('weevil_phyloseq.rds')
physeq <- readRDS('weevil_phyloseq.rds')
print(physeq)
plot_bar(physeq, fill='Phylum')
plot_bar(physeq.rarefied, fill='Phylum')
plot_bar(physeq.rarefied, x='Weevils', fill='Phylum')
plot_bar(physeq.rarefied, fill='Phylum')
physeq.rarefied <- readRDS('weevil_phyloseq_rarefied.rds')
print(physeq.rarefied)
plot_bar(physeq.rarefied, fill='Phylum')
plot_bar(physeq.rarefied, x='Weevils, fill='Phylum')
plot_bar(physeq.rarefied, x='', fill='Phylum')
