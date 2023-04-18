## This script goes through making a taxa bar plot and then running
## DESeq2 to find differentially abundant ASVs


# for help installing phyloseq, see this website
# https://bioconductor.org/packages/release/bioc/html/phyloseq.html

# to install phyloseq:
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("phyloseq")

library(qiime2R)
library(phyloseq)
library(zoo)
library(tidyverse)
devtools::install_github("r-lib/conflicted")

##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
##############################################

setwd('~/Desktop/output-all')
list.files()

if(!dir.exists("output"))
  dir.create("output")

metadata<-read.delim("~/Desktop/output-all/bbmetadata.txt")
str(metadata)
colnames(metadata)[3] = "culture"
colnames(metadata)[4] = "beater"
colnames(metadata)[5] = "time"
levels(metadata$beater) 

metadata$beater.ord = factor(metadata$beater, c("No Beating", "FastPrep", "BioSpec", "Positive Control"))
levels(metadata$culture)
metadata$culture.ord = factor(metadata$culture, c("Even", "High-gram-positive", "Positive Control"))
levels(metadata$time)
metadata$time.ord = factor(metadata$time, c("0 min", "10 min", "40 min"))
levels(metadata$beater.ord)
levels(metadata$culture.ord)
levels(metadata$time.ord)
#colnames(metadata)[8] <- "reported.antibiotic.usage"
#colnames(metadata)[9] <- "days.since.experiment.start"
str(metadata)

row.names(metadata) <- metadata[ ,1]

#metadata <- metadata[,-1]
row.names(metadata)

##Qiime2r method of reading in the taxonomy files
taxonomy<-read_qza("classifier-training/taxonomy.qza") 
head(taxonomy$data)

tax.clean<-parse_taxonomy(taxonomy$data)
head(tax.clean)

#All this is OK except that in future use of the taxonomy table, 
#these ASVs will be ignored because they are not classified. Why 
#are ASVs not classified? Its because there is not a close enough 
#match in the database. Just because there is not a good match in 
#the database does not mean they don’t exist, so I wanted to make 
#sure this data was not lost. So in my new code, from lines 200 – 224 
#I make it so that ASVs that are unclassified at any level are 
#classified as the lowest taxonomic level for which there is a 
#classification.
#Next, all these `NA` classifications with the last level that was 
#classified

tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

weight ="bold"
size = "size"


#################################################################
##Taxa barplot
#################################################################

physeq <- qza_to_phyloseq(
  features="core-metrics-results/rarefied_table.qza",
  tree="rooted-tree.qza",
  taxonomy = "classifier-training/taxonomy.qza",
  metadata = "bbmetadata.txt"
)



#First get the OTU table from physeq
physeq_otu_table <- data.frame(otu_table(physeq), check.names = F)

tax.clean = tax.clean[row.names(tax.clean) %in% rownames(physeq_otu_table),]
metadata.filtered = metadata[row.names(metadata) %in% colnames(physeq_otu_table),]

#Assign as variables to be feed into phyloseq
OTU.physeq = otu_table(as.matrix(physeq_otu_table), taxa_are_rows=TRUE)

#our edited and formatted taxonomy table from the top of this script
tax.physeq = tax_table(as.matrix(tax.clean))    
meta.physeq = sample_data(metadata.filtered)

#We then merge these into an object of class phyloseq.

physeq_bar_plot = phyloseq(OTU.physeq, tax.physeq, meta.physeq)



# Set colors for plotting
my_colors <- hcl(seq(15, 275, length.out = length(unique(physeq_meta_filtered[[ml]]))),
                 100, 65)
my_colors <- c( 'blue','#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)

#If you want different taxonomic level, find and replace the taxonomic level listed here
my_level <- c("Species")
my_column <- "beater" #this is the metadata column that we will use in the taxa barplot
my_column2 <- "culture"
my_column3 <- "time"

rm(taxa.summary)

abund_filter <- 0.01  # Our abundance threshold
ml = "Species" 
##====================
#taxon_abundance <- physeq_bar_plot %>%
#select(Taxon, Abundance, Time, Culture)


for(ml in my_level){
  print(ml)
  taxa.summary <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column),get(my_column2),get(my_column3), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  
colnames(taxa.summary)[2] <- my_column2
  colnames(taxa.summary)[3] <- my_column3
  colnames(taxa.summary)[4] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  
  physeq_meta_filtered$beater.ord = factor(physeq_meta_filtered$beater, c("No Beating", "BioSpec", "FastPrep", "Positive Control"))
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = get(my_column), y = Abundance.average, fill = get(ml))) + 
    facet_grid(time~culture) +
    geom_bar(stat = "identity") +
    #scale_fill_brewer(palette = "Paired") +
    scale_fill_manual(values = my_colors) +
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    theme_classic() +
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = F, keywidth = 1, keyheight = 1, ncol = 1)) +
    theme(legend.text=element_text(size=10)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=10)) +
    theme(legend.title = element_text(face=weight)) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 sample")) +
    labs(fill = "Diversity Index") +
    theme(axis.text.x = element_text(face = weight, size = 10),
       axis.text.y = element_text(face = weight, size = 10),
       axis.title = element_text(face = weight, size = 14),
       legend.text = element_text(face = weight),
       legend.title = element_text(face = weight, size = 14),
       plot.title = element_text(face = "bold", hjust = 0.5),
       strip.text = element_text(face = weight, size = '12')) +
        theme(panel.grid.major = element_line(size = 0.1, linetype = 'dotted', colour = "black"), 
        panel.grid.minor = element_line(size = 0.1, linetype = 'dotted', colour = "black")) # Increase the grid density
  ggsave(paste0("output/", ml, "Taxa-bar-plot", my_column2, ".png"), height = 5, width = 4, dpi=1000)
}


                    