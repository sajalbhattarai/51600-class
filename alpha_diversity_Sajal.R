
setwd("/Users/sajalbhattarai/Desktop/output-all/core-metrics-results")
list.files()


if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20

library(tidyverse)
library(qiime2R)
library(ggpubr)
#devtools::install_github("r-lib/conflicted")
#library(conflicted)
#library(dplyr)
meta <-read.delim("bbmetadata.txt")
str(meta)
meta$culture <- factor(meta$culture, levels = c("Even", "High-gram-positive", "Positive Control"))
meta$beater <- factor(meta$beater, levels = c("No Beating", "FastPrep", "BioSpec", "Positive Control"))
meta$time <- factor (meta$time, levels =c ("0 min", "10 min", "40 min"))

evenness <- read_qza("evenness_vector.qza")
evenness <- evenness$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

observed_features <- read_qza("observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

shannon <- read_qza("shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

faith_pd <- read_qza("faith_pd_vector.qza")
faith_pd <- faith_pd$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\
faith_pd <- faith_pd %>% dplyr::select(V1, V2)
colnames(faith_pd) <- c("SampleID","faith_pd")

###Alpha Diversity tables
# These tables will be merged for convenience and added to the 
# metadata table as the original tutorial was organized.


alpha_diversity <- merge(x=faith_pd, y=evenness, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, observed_features, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, shannon, by.x = "SampleID", by.y = "SampleID")
meta = merge(meta, alpha_diversity, by.x = "SampleID", by.y = "SampleID")
str(meta)

#Plots

#set the font details
font <- "Times New Roman"
size <- 12
weight <- "bold"
#rm(dpi)
par(bg="white", font.main=2, font.lab=2, cex.main=1, cex.lab=1)
hist(meta$shannon_entropy, main="Shannon diversity", xlab="", breaks=10)
hist(meta$faith_pd, main="Faith phylogenetic diversity", xlab="", breaks=10)
hist(meta$pielou_e, main="Evenness", xlab="", breaks=10)
hist(as.numeric(meta$observed_features), main="Observed Features", xlab="", breaks=10)

#Plots the qq-plot for residuals
theme(text = element_text(size = size, face = weight))
ggqqplot(meta$shannon_entropy, title = "Shannon-QQ plot") + theme(text = element_text(size = size, face = weight)) + theme (plot.title = element_text(size=size, hjust=0.5))
ggqqplot(meta$faith_pd, title = "Faith PD-QQ plot") + theme(text = element_text(size = size, face = weight)) + theme (plot.title = element_text(size=size, hjust=0.5))
ggqqplot(meta$pielou_e, title = "Evenness-QQ plot") + theme(text = element_text(size = size, face = weight)) + theme (plot.title = element_text(size=size, hjust=0.5))
ggqqplot(meta$observed_features, title = "Observed Features-QQ plot") + theme(text = element_text(size = size, face = weight)) + theme (plot.title = element_text(size=size, hjust=0.5))


shapiro.test(meta$shannon)
shapiro.test(meta$faith_pd)
shapiro.test(meta$pielou_e)
shapiro.test(meta$observed_features)
#=================

#===============
#EVENNESS
aov.evenness.beater = aov(pielou_evenness ~ beater, data=meta)
summary(aov.evenness.beater)

TukeyHSD(aov.evenness.beater)

levels(meta$beater)

#Plot
boxplot(pielou_evenness ~ beater, data=meta, ylab="Peilou evenness")
ggplot(meta, aes(beater, pielou_evenness)) + 
  geom_boxplot() + 
  facet_grid(cols = vars(culture), rows = vars(time)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = weight, size ="8"),
        axis.title = element_text(face = weight, size = "10"),
        strip.text = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "10"),
        legend.title = element_text(face = weight, size = "10"), axis.text.y = element_text(face=weight, size = "8"))
#========================

#====================

#SHANNON_ENTROPY shannon_entropy
aov.shannon.beater = aov(shannon_entropy ~ beater, data=meta)
summary(aov.shannon.beater)

TukeyHSD(aov.shannon.beater)

levels(meta$beater)

#aov.shannon = aov(shannon_entropy ~ time + culture + beater, data = meta)
#summary(aov.shannon)
#TukeyHSD(aov.shannon, "beater")

#Plot
boxplot(shannon_entropy ~ beater, data=meta, ylab="Shannon")
ggplot(meta, aes(beater, shannon_entropy)) + 
  geom_boxplot() + 
  facet_grid(cols = vars(culture), rows = vars(time)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = weight, size ="8"),
        axis.title = element_text(face = weight, size = "10"),
        strip.text = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "10"),
        legend.title = element_text(face = weight, size = "10"), axis.text.y = element_text(face=weight, size = "8"))
#===================
#OBSERVED FEATURES
aov.observed_features.beater = aov(observed_features ~ beater, data=meta)
summary(aov.observed_features.beater)

TukeyHSD(aov.observed_features.beater)

levels(meta$beater)

#Plot
boxplot(observed_features ~ beater, data=meta, ylab="Observed features")
ggplot(meta, aes(beater, observed_features)) + 
  geom_boxplot() + 
  facet_grid(cols = vars(culture), rows = vars(time)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = weight, size ="8"),
        axis.title = element_text(face = weight, size = "10"),
        strip.text = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "10"),
        legend.title = element_text(face = weight, size = "10"), axis.text.y = element_text(face=weight, size = "8"))
#===================
#FAITH_PD
aov.faith_pd.beater = aov(faith_pd ~ beater, data=meta)
summary(aov.faith_pd.beater)

TukeyHSD(aov.faith_pd.beater)

levels(meta$beater)
boxplot(faith_pd ~ beater, data=meta, ylab="Faith PD")
ggplot(meta, aes(beater, faith_pd)) + 
  geom_boxplot() + 
  facet_grid(cols = vars(culture), rows = vars(time)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = weight, size ="8"),
        axis.title = element_text(face = weight, size = "10"),
        strip.text = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "10"),
        legend.title = element_text(face = weight, size = "10"), axis.text.y = element_text(face=weight, size = "8"))


#========================
# ggsave("output/evenness_boxplot.png", evenness, height = 3, width = 3)
# =======================

# Now, the above graph is kind of not correct. Our test and our graphic do not exactly match. ANOVA and Tukey are tests based on the mean, but the boxplot plots the median. Its not wrong, its just not the best method. Unfortunately plotting the average and standard deviation is a little complicated.
#Mean Plots
##=========================
#Evenness
evenness_summary <- meta %>% 
  group_by(beater, time, culture) %>% # add time and culture as grouping variables
  summarise(mean_evenness = mean(pielou_evenness),  
            sd_evenness = sd(pielou_evenness), 
            n_evenness = n(),  
            se_evenness = sd(pielou_evenness)/sqrt(n())) 

ggplot(evenness_summary, aes(beater, mean_evenness, fill = beater)) + 
  geom_boxplot() +
  geom_errorbar(aes(ymin = mean_evenness - se_evenness, ymax = mean_evenness + se_evenness), size = 0.1, linewidth = 0.25) + 
  facet_grid(rows = vars(time), cols = vars(culture)) + # add time and culture axes
  theme_q2r() +
  labs(y="Evenness  ± s.e.", x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = weight, size ="8"),
        axis.title = element_text(face = weight, size = "10"),
        strip.text = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "10"),
        legend.title = element_text(face = weight, size = "10"), axis.text.y = element_text(face=weight, size = "8"))
  
#theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3, face =weight)) +
  theme(axis.text.y = element_text(size = 3, face = weight, hjust = 1)) + # adjust font size of y-axis labels
  theme(legend.title = element_text(size = 3, face = weight), legend.text = element_text(size = 3, face = weight)) + # adjust font size of legend text) +
  labs(y="Pielou's evenness  ± s.e.", x = "") +
  theme(plot.title = element_text(size=1, face =weight)) + # adjust font size of title
  theme(strip.text = element_text(size=3, face =weight)) # adjust font size of facet labels
  #theme(axis.title.y = element_text(size =3, face = weight))
ggsave("evenness_se.png", evenness_se, height = 2.5, width = 3, dpi = 800)

##
#SHANNON
shannon_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(beater, time, culture) %>%   # the grouping variable
  summarise(mean_shannon = mean(shannon_entropy),  # calculates the mean of each group
            sd_shannon = sd(shannon_entropy), # calculates the standard deviation of each group
            n_shannon = n(),  # calculates the sample size per group
            se_shannon = sd(shannon_entropy)/sqrt(n())) # calculates the standard error of each group
# We can now make a bar plot of means vs body site, with standard 
# deviations or standard errors as the error bar. The following code 
# uses the standard deviations.

ggplot(shannon_summary, aes(beater, mean_shannon, fill = beater)) + 
  geom_boxplot() +
  geom_errorbar(aes(ymin = mean_shannon - se_shannon, ymax = mean_shannon + se_shannon), size = 0.1, linewidth = 0.25) + 
  facet_grid(rows = vars(time), cols = vars(culture)) + # add time and culture axes
  theme_q2r() +
  labs(y="Shannon  ± s.e.", x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = weight, size ="8"),
        axis.title = element_text(face = weight, size = "10"),
        strip.text = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "10"),
        legend.title = element_text(face = weight, size = "10"), axis.text.y = element_text(face=weight, size = "8"))

#theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3, face =weight)) +
  theme(axis.text.y = element_text(size = 3, face = weight, hjust = 1)) + # adjust font size of y-axis labels
  theme(legend.title = element_text(size = 3, face = weight), legend.text = element_text(size = 3, face = weight)) + # adjust font size of legend text) +
  labs(y="Shannon  ± s.e.", x = "") +
  theme(plot.title = element_text(size=1, face =weight)) + # adjust font size of title
  theme(strip.text = element_text(size=3, face =weight)) # adjust font size of facet labels
ggsave("shannon_se.png", shannon_se, height = 2.5, width = 3, dpi = 800)
##=========================
#observed_features
observed_features_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(beater, time, culture) %>%   # the grouping variable
  summarise(mean_observed_features = mean(observed_features),  # calculates the mean of each group
            sd_observed_features = sd(observed_features), # calculates the standard deviation of each group
            n_observed_features = n(),  # calculates the sample size per group
            se_observed_features = sd(observed_features)/sqrt(n()))
# calculates the standard error of each group
# We can now make a bar plot of means vs body site, with standard 
# deviations or standard errors as the error bar. The following code 
# uses the standard deviations.

ggplot(observed_features_summary, aes(beater, mean_observed_features, fill = beater)) + 
  geom_boxplot() +
  geom_errorbar(aes(ymin = mean_observed_features - se_observed_features, ymax = mean_observed_features + se_observed_features), size = 0.1, linewidth = 0.25) + 
  facet_grid(rows = vars(time), cols = vars(culture)) + # add time and culture axes
  theme_q2r() +
  labs(y="Observed features  ± s.e.", x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = weight, size ="8"),
        axis.title = element_text(face = weight, size = "10"),
        strip.text = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "10"),
        legend.title = element_text(face = weight, size = "10"), axis.text.y = element_text(face=weight, size = "8"))
  
#theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3, face =weight)) +
  theme(axis.text.y = element_text(size = 3, face = weight, hjust = 1)) + # adjust font size of y-axis labels
  theme(legend.title = element_text(size = 3, face = weight), legend.text = element_text(size = 3, face = weight)) + # adjust font size of legend text) +
  labs(y="observed_features  ± s.e.", x = "") +
  theme(plot.title = element_text(size=1, face =weight)) + # adjust font size of title
  theme(strip.text = element_text(size=3, face =weight)) # adjust font size of facet labels
ggsave("observed_features_se.png", observed_features_se, height = 2.5, width = 3, dpi = 800)
##===================
#faith_pd
faith_pd_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(beater, time, culture) %>%   # the grouping variable
  summarise(mean_faith_pd = mean(faith_pd),  # calculates the mean of each group
            sd_faith_pd = sd(faith_pd), # calculates the standard deviation of each group
            n_faith_pd = n(),  # calculates the sample size per group
            se_faith_pd = sd(faith_pd)/sqrt(n()))
# calculates the standard error of each group
# We can now make a bar plot of means vs body site, with standard 
# deviations or standard errors as the error bar. The following code 
# uses the standard deviations.

ggplot(faith_pd_summary, aes(beater, mean_faith_pd, fill = beater)) + 
  geom_boxplot() +
  geom_errorbar(aes(ymin = mean_faith_pd - se_faith_pd, ymax = mean_faith_pd + se_faith_pd), size = 0.1, linewidth = 0.25) + 
  facet_grid(rows = vars(time), cols = vars(culture)) + # add time and culture axes
  theme_q2r() +
  labs(y="Faith PD  ± s.e.", x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = weight, size ="8"),
        axis.title = element_text(face = weight, size = "10"),
        strip.text = element_text(face = weight, size = "10"),
        legend.text = element_text(face = weight, size = "10"),
        legend.title = element_text(face = weight, size = "10"), axis.text.y = element_text(face=weight, size = "8"))
  
  
#theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3, face =weight)) +
  theme(axis.text.y = element_text(size = 3, face = weight, hjust = 1)) + # adjust font size of y-axis labels
  theme(legend.title = element_text(size = 3, face = weight), legend.text = element_text(size = 3, face = weight)) + # adjust font size of legend text) +
  labs(y="faith_pd  ± s.e.", x = "") +
  theme(plot.title = element_text(size=1, face =weight)) + # adjust font size of title
  theme(strip.text = element_text(size=3, face =weight)) # adjust font size of facet labels
ggsave("faith_pd_se.png", faith_pd_se, height = 2.5, width = 3, dpi = 800)

####=============
## **Non-normally distributed metrics**

# We will use **Faith's phylogenetic diversity** here. Since body site 
# is categorical, we use Kruskal-Wallis (non-parametric equivalent of 
# ANOVA). If we have only two levels, we would run Wilcoxon rank sum 
# test (non-parametric equivalent of t-test)

kruskal.test(faith_pd ~ beater, data=meta)

# We can test pairwise within the age groups with Wilcoxon Rank Sum 
# Tests. This test has a slightly different syntax than our other tests

pairwise.wilcox.test(meta$faith_pd, meta$beater, p.adjust.method="BH")

# Like evenness, we see that pd also increases with age.

#Plot
boxplot(faith_pd ~ beater, data=meta, ylab="Faith phylogenetic diversity")

ggplot(meta, aes(x = beater, y = faith_pd)) +
  geom_boxplot(aes(color=beater)) +
  facet_grid(rows = vars(time), cols = vars(culture)) +
  ylab("Faith phylogenetic diversity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# or with ggplot2

ggplot(meta, aes(beater, faith_pd)) + 
  geom_boxplot(aes(color = beater)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith Phylogenetic Diversity", x = "") 


##Continuous variables
# For continuous variables, we use general linear models, specifying 
# the distribution that best fits our data.

# **Normally distributed metrics**

# Since days.since.experiment.start is a continuous variable, we run a 
# general linear model. We will again use evenness as our roughly normal 
# metric. The default of `glm` and `lm` is the normal distribution so we 
# don't have to specify anything.

# Does days.since.experiment.start impact evenness of the microbiota?

glm.evenness.time = glm(pielou_evenness ~ time, data=meta)
summary(glm.evenness.time)

#The output let's us know that the intercept of our model is significantly different from 0 but our slope (*e.g.* our variable of interest) is not. This makes sense when we look at the data.

plot(pielou_evenness ~ time, data=meta)
#Add the glm best fit line
plot(pielou_evenness ~ time, data=meta) + abline(glm.evenness.time)

# **Non-normally distributed metrics**

# We will again use a *general linear model* for our non-normally 
# distributed metric Faith_pd. However, this time, we change the 
# distribution from normal to something that fits the data better. 

# But which distribution should we choose? In statistics, there is no 
# one "best" model. There are only good and better models. We will use 
# the plot() function to compare two models and pick the better one.

# First, the Gaussian (normal) distribution, which we already know is a bad fit.

gaussian.faith.time = glm(faith_pd ~ time, data=meta, family="gaussian")
plot(gaussian.faith.time, which=c(1,2))

# Quasipoisson (log) distribution
qp.faith.time = glm(faith_pd ~ time, data=meta, family="quasipoisson")
plot(qp.faith.time, which=c(1,2))

# What we're looking for is no pattern in the Residuals vs. Fitted graph 
# ("stars in the sky"), which shows that we picked a good distribution 
# family to fit our data. We also want our residuals to be normally 
# distributed, which is shown by most/all of the points falling on the 
# line in the Normal Q-Q plot.

# While it's still not perfect, the quasipoisson fits much better. 
# In the residuals vs fitted graph, the y axis is from -2 to 4  whereas 
# the axis with gaussian was from -5 to 10. So, we will use quasipoisson 
# and see that ADG does not to correlate to Chao richness.
summary(qp.faith.time)

# Plotting this we see that, indeed, there is a trend toward correlation between Faith_pd and time.

#Plot
plot(log(faith_pd) ~ time, data=meta, ylab="ln(Faith Phylo. Diversity)")
plot(log(faith_pd) ~ time, data=meta, ylab="ln(Faith Phylo. Diversity)") + abline(qp.faith.time)

