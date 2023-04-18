
#Load the packages
install.packages("tidyverse")
install.packages("vegan")
install.packages("devtools")
library(devtools)
devtools::install_github("jbisanz/qiime2R")


library(tidyverse)
library(vegan)
library(qiime2R)
###Set your working directory
#
weight = "bold"
size = "size"
setwd("/Users/sajalbhattarai/Desktop/output-all/core-metrics-results")

list.files()

if(!dir.exists("beta-diversity-output"))
  dir.create("beta-diversity-output")

#Now the qiime2R method
metadata <- read.delim("~/Desktop/output-all/core-metrics-results/bbmetadata.txt") 
str(metadata)
# levels(metadata$`body-site`)
# colnames(metadata)[3] <- "beater"
# colnames(metadata)[8] <- "reported.antibiotic.usage"
# colnames(metadata)[9] <- "days.since.experiment.start"
# str(metadata)

row.names(metadata) <- metadata[,1] 
row.names(metadata) <- metadata$SampleID
#metadata <- metadata[,-1]
row.names(metadata)

bc_PCoA<-read_qza("bray_curtis_pcoa_results.qza")
wUF <- read_qza("weighted_unifrac_pcoa_results.qza")
uwUF <- read_qza("unweighted_unifrac_pcoa_results.qza")
jac <- read_qza("jaccard_pcoa_results.qza")
body_colors <- c("orange", "Blue","Black", "Red")
body_colors2 <- c( "Blue","Black", "Red", "orange")

##==============

bc_meta <- bc_PCoA$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

wUF_meta <- wUF$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

uwUF_meta <- uwUF$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

jac_meta <- jac$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))
##===============

# Now we are going to make an ordination plot
# ggplot(bc_meta, aes(x=PC1, y=PC2, color=beater)) +
#   geom_point() + #alpha controls transparency and helps when points are overlapping
#   theme_q2r() +
#   xlab("PC1 (32.27%)") +
#   ylab("PC2 (22.28%)") +
#   scale_color_manual(values=c("Blue", "Black", "Green", "Gray"), name = "body-site")

# Now we are going to make our code a little more re-usable
my_column <- "beater"
my_column2 <- "culture"
my_column3 <- "time"
#my_column <- "DietTreatment"
####=================
###BRAY_CURTIS
ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  facet_grid(culture~time) +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column) +
  ggtitle("Bray-Curtis-PCOA-results") +
  theme(axis.text.x = element_text(face = weight, size = "10"),
        axis.text.y = element_text(face = weight, size = "10"),
        axis.title = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "8"),
        legend.title = element_text(face = weight, size = "8"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = "14"),
        strip.text = element_text(face = weight, size = "10"))
ggsave(paste0("beta-diversity-output/BC-basic_", my_column,".tiff"), height=3, width=5, device="tiff", dpi = 800) # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "beater"
weight <- "bold"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= culture), size = 1.5) +
  scale_shape_discrete(name = "Culture shape") +
  #scale_color_manual(values = body_colors, name = "Beater color")+
  #scale_color_discrete(name = "Beater color") + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 2) +
  theme_q2r() +
  facet_grid(~time) +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors2, name = "Beater color") +
  ggtitle("Centroid-Bray-Curtis-PCOA-results") +
  theme(axis.text.x = element_text(face = weight, size = "10"),
        axis.text.y = element_text(face = weight, size = "10"),
        axis.title = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "8"),
        legend.title = element_text(face = weight, size = "8"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = "14"),
        strip.text = element_text(face = weight, size = "10"))

#######===============
###weighted Unifrac
ggplot(wUF_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  facet_grid(culture~time) +
  xlab(paste0("PC1 (", round(100*wUF$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*wUF$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column) +
  ggtitle("Weighted Unifrac-PCOA-results") +
  theme(axis.text.x = element_text(face = weight, size = "10"),
        axis.text.y = element_text(face = weight, size = "10"),
        axis.title = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "8"),
        legend.title = element_text(face = weight, size = "8"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = "14"),
        strip.text = element_text(face = weight, size = "10"))
ggsave(paste0("beta-diversity-output/BC-basic_", my_column,".tiff"), height=3, width=5, device="tiff", dpi = 800) # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),wUF_meta,mean)
colnames(centroids)[1] <- "beater"
weight <- "bold"

ggplot(wUF_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= culture), size = 1.5) +
  scale_shape_discrete(name = "Culture shape") +
  #scale_color_manual(values = body_colors, name = "Beater color")+
  #scale_color_discrete(name = "Beater color") + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 2) +
  theme_q2r() +
  facet_grid(~time) +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*wUF$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*wUF$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors2, name = "Beater color") +
  ggtitle("Centroid-Weighted Unifrac-PCOA-results") +
  theme(axis.text.x = element_text(face = weight, size = "10"),
        axis.text.y = element_text(face = weight, size = "10"),
        axis.title = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "8"),
        legend.title = element_text(face = weight, size = "8"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = "14"),
        strip.text = element_text(face = weight, size = "10"))
#==================
#Unweighted Unifrac
ggplot(uwUF_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  facet_grid(culture~time) +
  xlab(paste0("PC1 (", round(100*uwUF$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*uwUF$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column) +
  ggtitle("Unweighted Unifrac-PCOA-results") +
  theme(axis.text.x = element_text(face = weight, size = "10"),
        axis.text.y = element_text(face = weight, size = "10"),
        axis.title = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "8"),
        legend.title = element_text(face = weight, size = "8"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = "14"),
        strip.text = element_text(face = weight, size = "10"))
ggsave(paste0("beta-diversity-output/BC-basic_", my_column,".tiff"), height=3, width=5, device="tiff", dpi = 800) # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),uwUF_meta,mean)
colnames(centroids)[1] <- "beater"
weight <- "bold"

ggplot(uwUF_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= culture), size = 1.5) +
  scale_shape_discrete(name = "Culture shape") +
  #scale_color_manual(values = body_colors, name = "Beater color")+
  #scale_color_discrete(name = "Beater color") + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 2) +
  theme_q2r() +
  facet_grid(~time) +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*uwUF$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*uwUF$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors2, name = "Beater color") +
  ggtitle("Centroid-Unweighted Unifrac-PCOA-results") +
  theme(axis.text.x = element_text(face = weight, size = "10"),
        axis.text.y = element_text(face = weight, size = "10"),
        axis.title = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "8"),
        legend.title = element_text(face = weight, size = "8"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = "14"),
        strip.text = element_text(face = weight, size = "10"))
###
####
####
#Jaccard
ggplot(jac_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  facet_grid(culture~time) +
  xlab(paste0("PC1 (", round(100*jac$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*jac$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column) +
  ggtitle("Jaccard-PCOA-results") +
  theme(axis.text.x = element_text(face = weight, size = "10"),
        axis.text.y = element_text(face = weight, size = "10"),
        axis.title = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "8"),
        legend.title = element_text(face = weight, size = "8"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = "14"),
        strip.text = element_text(face = weight, size = "10"))
ggsave(paste0("beta-diversity-output/BC-basic_", my_column,".tiff"), height=3, width=5, device="tiff", dpi = 800) # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),jac_meta,mean)
colnames(centroids)[1] <- "beater"
weight <- "bold"

ggplot(jac_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= culture), size = 1.5) +
  scale_shape_discrete(name = "Culture shape") +
  #scale_color_manual(values = body_colors, name = "Beater color")+
  #scale_color_discrete(name = "Beater color") + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 2) +
  theme_q2r() +
  facet_grid(~time) +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*jac$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*jac$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors2, name = "Beater color") +
  ggtitle("Centroid-Jaccard-PCOA-results") +
  theme(axis.text.x = element_text(face = weight, size = "10"),
        axis.text.y = element_text(face = weight, size = "10"),
        axis.title = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "8"),
        legend.title = element_text(face = weight, size = "8"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = "14"),
        strip.text = element_text(face = weight, size = "10"))
########=========
#PERMANOOVA TESTS
##BRAY_CURTIS
bc_dist_mat<-read_qza("bray_curtis_distance_matrix.qza")
bc_dm <- as.matrix(bc_dist_mat$data) 
#length(rownames(bc_dm))
#length(sample_ids)
#setdiff(rownames(bc_dm), sample_ids)
#setdiff(sample_ids, rownames(bc_dm))
#metadata2 <- subset(metadata, SampleID %in% rownames(bc_dm))
rownames(bc_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_bc_sub <- metadata[match(rownames(bc_dm),metadata$SampleID),]
rownames(bc_dm) == metadata_bc_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(bc_dm ~ beater, data = metadata_bc_sub)

write.table(PERMANOVA_out,"bc_Beater_Adonis_overall.csv",sep=",", row.names = TRUE) 

######################################################################################
##  Pairwise adonis function
##  we can also performe a pairwise comparison with the function 
##  Pairwise Adonis funtion by edro Martinez Arbizu & Sylvain Monteux
##  https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
#######################################################################################

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 

beater_Pair <- pairwise.adonis2(bc_dm ~ beater, data = metadata_bc_sub)
write.table(beater_Pair,"bc_beater_Adonis_pairwise.csv",sep=",", row.names = TRUE) 

####============
#WEIGHTED_UNIFRAC
wuF_dist_mat<-read_qza("weighted_unifrac_distance_matrix.qza")
wuF_dm <- as.matrix(wuF_dist_mat$data) 
#length(rownames(wuF_dm))
#length(sample_ids)
#setdiff(rownames(wuF_dm), sample_ids)
#setdiff(sample_ids, rownames(wuF_dm))
#metadata2 <- subset(metadata, SampleID %in% rownames(wuF_dm))
rownames(wuF_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_wuF_sub <- metadata[match(rownames(wuF_dm),metadata$SampleID),]
rownames(wuF_dm) == metadata_wuF_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(wuF_dm ~ beater, data = metadata_wuF_sub)

write.table(PERMANOVA_out,"wuF_Beater_Adonis_overall.csv",sep=",", row.names = TRUE) 

######################################################################################
##  Pairwise adonis function
##  we can also performe a pairwise comparison with the function 
##  Pairwise Adonis funtion by edro Martinez Arbizu & Sylvain Monteux
##  https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
#######################################################################################

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 

beater_Pair <- pairwise.adonis2(wuF_dm ~ beater, data = metadata_wuF_sub)
write.table(beater_Pair,"wuF_beater_Adonis_pairwise.csv",sep=",", row.names= TRUE)

########
##UNWEIGHTED_UNIFRAC
uwUF_dist_mat<-read_qza("unweighted_unifrac_distance_matrix.qza")
uwUF_dm <- as.matrix(uwUF_dist_mat$data) 
#length(rownames(uwUF_dm))
#length(sample_ids)
#setdiff(rownames(uwUF_dm), sample_ids)
#setdiff(sample_ids, rownames(uwUF_dm))
#metadata2 <- subset(metadata, SampleID %in% rownames(uwUF_dm))
rownames(uwUF_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_uwUF_sub <- metadata[match(rownames(uwUF_dm),metadata$SampleID),]
rownames(uwUF_dm) == metadata_uwUF_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(uwUF_dm ~ beater, data = metadata_uwUF_sub)

write.table(PERMANOVA_out,"uwUF_Beater_Adonis_overall.csv",sep=",", row.names = TRUE) 

######################################################################################
##  Pairwise adonis function
##  we can also performe a pairwise comparison with the function 
##  Pairwise Adonis funtion by edro Martinez Arbizu & Sylvain Monteux
##  https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
#######################################################################################

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 

beater_Pair <- pairwise.adonis2(uwUF_dm ~ beater, data = metadata_uwUF_sub)
write.table(beater_Pair,"uwUF_beater_Adonis_pairwise.csv",sep=",", row.names = TRUE) 

####=========
##JACCARD

jac_dist_mat<-read_qza("jaccard_distance_matrix.qza")
jac_dm <- as.matrix(jac_dist_mat$data) 
#length(rownames(jac_dm))
#length(sample_ids)
#setdiff(rownames(jac_dm), sample_ids)
#setdiff(sample_ids, rownames(jac_dm))
#metadata2 <- subset(metadata, SampleID %in% rownames(jac_dm))
rownames(jac_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_jac_sub <- metadata[match(rownames(jac_dm),metadata$SampleID),]
rownames(jac_dm) == metadata_jac_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(jac_dm ~ beater, data = metadata_jac_sub)

write.table(PERMANOVA_out,"jac_Beater_Adonis_overall.csv",sep=",", row.names = TRUE) 

######################################################################################
##  Pairwise adonis function
##  we can also performe a pairwise comparison with the function 
##  Pairwise Adonis funtion by edro Martinez Arbizu & Sylvain Monteux
##  https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
#######################################################################################

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 

beater_Pair <- pairwise.adonis2(jac_dm ~ beater, data = metadata_jac_sub)
write.table(beater_Pair,"jac_beater_Adonis_pairwise.csv",sep=",", row.names = TRUE) 
