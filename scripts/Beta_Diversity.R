#-------------------------#
####---Load Packages---####
#-------------------------#

library(pairwiseAdonis)
library(AICcPermanova)
library(lme4)
library(lsmeans)
library(lmerTest)
library(ggplot2)
library(ggpubr)
library(ggnewscale)
library(gridExtra)
library(phyloseq)
library(microViz) #ps_filter function 

#---------------------#
####---Load Data---####
#---------------------#

data_fd <- "results/3 - Create Phyloseq/"

results_fd <- "results/6 - Beta-diversity/"
if(! dir.exists(results_fd)) { dir.create(results_fd, recursive = TRUE) }


mybiom2 <- readRDS(paste0(data_fd, "mybiom2.rds"))

source("scripts/Functions.R")

#-------------------------#
####---Set variables---####
#-------------------------#

# Color palettes
cohort_cols <- c("#2C7800", "#b54e5d")
tissue_cols2 <- c("#ca5686", "#ca7040")

# Unique groups
uni_tiss <- c("Hemolymph", "Hepatopancreas")
uni_days <- day_levels <- c("0", "7", "14", "21", "28", "35", "37", "39")

# Distance methods
dist_meth <- c("bray", "wunifrac", "unifrac", "jaccard") 

# Set seed for reproducibility
my_seed <- 123  

#------------------------------------------#
####---Subset mybiom object by tissue---####
#------------------------------------------#

# Hemolymph
mybiom_hm <- ps_filter(mybiom2, Tissue == "Hemolymph")
sample_data(mybiom_hm)[, "Day"] <- as.factor(get_variable(sample_data(mybiom_hm))[, "Day"])

# Hepatopancreas
mybiom_hp <- ps_filter(mybiom2, Tissue == "Hepatopancreas")
sample_data(mybiom_hp)[, "Day"] <- as.factor(get_variable(sample_data(mybiom_hp))[, "Day"])

#---------------------------#
####---Distance models---####
#---------------------------#

# Rarefy data
mybiom.r <- rarefy_even_depth(mybiom2, rngseed = my_seed, verbose = FALSE)
mybiom.hm.r <- rarefy_even_depth(mybiom_hm, rngseed = my_seed, verbose = FALSE)
mybiom.hp.r <- rarefy_even_depth(mybiom_hp, rngseed = my_seed, verbose = FALSE)

# Calculate distances 
## Empty lists for saving results
### Hemolymph
hm.bray <- list() 
hm.wuni <- list()

### Hepatopancreas
hp.bray <- list()
hp.wuni <- list()

# Permanova & Beta-dispersion tets
hm.bray <- beta_distance("Hemolymph", "bray", mybiom.hm.r)
hm.wuni <- beta_distance("Hemolymph", "wunifrac", mybiom.hm.r)

hp.bray <- beta_distance("Hepatopancreas", "bray", mybiom.hp.r)
hp.wuni <- beta_distance("Hepatopancreas", "wunifrac", mybiom.hp.r)

## Create figure plot
hm.bray.plot <- hm.bray$Plot + labs(tag = "A")
hm.wuni.plot <- hm.wuni$Plot + labs(tag = "B")
hp.bray.plot <- hp.bray$Plot + labs(tag = "C")
hp.wuni.plot <- hp.wuni$Plot + labs(tag = "D")

pdf(paste0(results_fd, "SuppFigure4.pdf"), width = 20, height = 20)

ggarrange(hm.bray.plot, hm.wuni.plot, hp.bray.plot, hp.wuni.plot, nrow = 2, ncol = 2)

dev.off()

# Modelling
## Empty lists for saving results
### Hemolymph
hm.bray.day.model <- list() 
hm.bray.cohort.model <- list()

hm.wuni.day.model <- list()
hm.wuni.cohort.model <- list()

### Hepatopancreas
hp.bray.day.model <- list()
hp.bray.cohort.model <- list()

hp.wuni.day.model <- list()
hp.wuni.cohort.model <- list()

## Run modelling
### Hemolymph
hm.bray.day.model <- beta_model("Hm", "bray", "Day", mybiom.hm.r)
hm.bray.cohort.model <- beta_model("Hm", "bray", "Cohort", mybiom.hm.r)

hm.wuni.day.model <- beta_model("Hm", "wunifrac", "Day", mybiom.hm.r)
hm.wuni.cohort.model <- beta_model("Hm", "wunifrac", "Cohort", mybiom.hm.r)

### Hepatopancreas
hp.bray.day.model <- beta_model("Hp", "bray", "Day", mybiom.hp.r)
hp.bray.cohort.model <- beta_model("Hp", "bray", "Cohort", mybiom.hp.r)

hp.wuni.day.model <- beta_model("Hp", "wunifrac", "Day", mybiom.hp.r)
hp.wuni.cohort.model <- beta_model("Hp", "wunifrac", "Cohort", mybiom.hp.r)

## Create figure plots
### Cohort
hm.b.c <- hm.bray.cohort.model$Plot + labs(tag = "A")
hm.w.c <- hm.wuni.cohort.model$Plot + labs(tag = "B")
hp.b.c <- hp.bray.cohort.model$Plot + labs(tag = "C")
hp.w.c <- hp.wuni.cohort.model$Plot + labs(tag = "D")

pdf(paste0(results_fd, "Figure4.pdf"), width = 10, height = 10)

ggarrange(hm.b.c, hm.w.c, hp.b.c, hp.w.c, nrow = 2, ncol = 2, common.legend = TRUE, legend = "right")

dev.off()

### Day (Supplementary)
hm.b.d <- hm.bray.day.model$Plot + labs(tag = "A")
hm.w.d <- hm.wuni.day.model$Plot + labs(tag = "B")
hp.b.d <- hp.bray.day.model$Plot + labs(tag = "C")
hp.w.d <- hp.wuni.day.model$Plot + labs(tag = "D")

pdf(paste0(results_fd, "Day_Supplementary.pdf"), width = 10, height = 10)

ggarrange(hm.b.d, hm.w.d, hp.b.d, hp.w.d, nrow = 2, ncol = 2, common.legend = TRUE, legend = "right")

dev.off()

#--------------------------------#
####---Distance Ordinations---####
#--------------------------------#

# Set seed for reproducibility <= ~~~~~~ FC: What is this for??? ~~~~~~ 
set.seed(my_seed)

# Create distance groups
## Hemolymph
dist_gr.hm <- ps_filter(mybiom.hm.r)
sample_data(dist_gr.hm)[, "Day"] <- as.factor(get_variable(sample_data(dist_gr.hm))[, "Day"]) 

## Hepatopancreas
dist_gr.hp <- ps_filter(mybiom.hp.r)
sample_data(dist_gr.hp)[, "Day"] <- as.factor(get_variable(sample_data(dist_gr.hp))[, "Day"]) 


#--------------------#
####---Plotting---####
#--------------------#

# Bray Curtis 
pdf(paste0(results_fd, "Bray_PCoA.pdf"), width = 10, height = 8)

## Hemolymph
bray.hm <- plotOrd(mybiom.hm.r, "bray", dist_gr.hm, color_var = "Day", shape_var = "Cohort", 
                   ellipse_var = "Day", plot_title = "Hemolymph - Bray Curtis") 

## Hepatopancreas
bray.hp <- plotOrd(mybiom.hp.r, "bray", dist_gr.hp, color_var = "Day", shape_var = "Cohort", 
                   ellipse_var = "Day", plot_title = "Hepatopancreas - Bray Curtis") 

ggarrange(bray.hm, bray.hp, nrow = 2, ncol = 1, common.legend = TRUE, legend = "bottom")

dev.off()

# UniFrac
pdf(paste0(results_fd, "UniFrac_PCoA.pdf"), width = 10, height = 8)

## Hemolymph
uni.hm <- plotOrd(mybiom.hm.r, "unifrac", dist_gr.hm, color_var = "Day", shape_var = "Cohort", 
                  ellipse_var = "Day", plot_title = "Hemolymph - UniFrac") 

## Hepatopancreas
uni.hp <- plotOrd(mybiom.hp.r, "unifrac", dist_gr.hp, color_var = "Day", shape_var = "Cohort", 
                  ellipse_var = "Day", plot_title = "Hepatopancreas - UniFrac") 

ggarrange(uni.hm, uni.hp, nrow = 2, ncol = 1, common.legend = TRUE, legend = "bottom")

dev.off()

# Weighted UniFrac
pdf(paste0(results_fd, "WUniFrac_PCoA.pdf"), width = 10, height = 8)

## Hemolymph
wuni.hm <- plotOrd(mybiom.hm.r, "wunifrac", dist_gr.hm, color_var = "Day", shape_var = "Cohort", 
                   ellipse_var = "Day", plot_title = "Hemolymph - W.UniFrac") 

## Hepatopancreas
wuni.hp <- plotOrd(mybiom.hp.r, "wunifrac", dist_gr.hp, color_var = "Day", shape_var = "Cohort", 
                   ellipse_var = "Day", plot_title = "Hepatopancreas - W.UniFrac") 

ggarrange(wuni.hm, wuni.hp, nrow = 2, ncol = 1, common.legend = TRUE, legend = "bottom")

dev.off()
