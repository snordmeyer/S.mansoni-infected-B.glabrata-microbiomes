#---------------------------#
####---Import Packages---####
#---------------------------#

library(microbiome) 
library(picante)
library(dplyr)
library(tidyr)
library(FSA)
library(lme4)
library(lsmeans)
library(ggpubr)

#-------------------#
####---Load Data---###
#--------------------#
data_fd <- "results/3 - Create Phyloseq/"

results_fd <- "results/5 - Alpha diversity/"

# Load objects
mybiom2 <- readRDS(paste0(data_fd, "mybiom2.rds"))
tree <- readRDS(paste0(data_fd, "tree.rds"))
asv <- readRDS(paste0(data_fd, "asv.rds"))
new_cols <- readRDS(paste0(data_fd, "new_cols.rds"))

source("scripts/Functions.R")

#-------------------------#
####---Set Variables---####
#-------------------------#

# Set levels
day_levels <- c("0", "7", "14", "21", "28", "35", "37", "39")
tissue_levels <- c("Hemolymph", "Hepatopancreas", "Water")

# Set colors
cohort_cols <- c("#2C7800", "#b54e5d")

#---------------------------#
####---Alpha Diversity---####
#---------------------------#

if(! dir.exists(results_fd)) { dir.create(results_fd, recursive = TRUE) }

# Create myalpha object
myalpha_index <- c("observed", "evenness_simpson", "PD")
myalpha <- microbiome::alpha(mybiom2, index = myalpha_index)
mypd    <- pd(t(asv), tree)
myalpha <- merge(myalpha, mypd, by="row.names")
myalpha <- merge(myalpha, new_cols[,1:11], by.x=1, by.y="row.names")

## Create myalpha subset containing only snail samples
myalpha.s <- myalpha %>% filter(Tissue != "Water")

# Set levels of Day variable 
myalpha$Day <- as.factor(myalpha$Day)
levels(myalpha$Day) <- day_levels

# Pull and save unique variables 
unique_tissues <- unique(myalpha$Tissue)
unique_CT <- unique(myalpha$Coh_Tis)
unique_CTD  <- unique(myalpha$Coh_Tis_Day)


unique_TD <- unique(myalpha$Tis_Day)
unique_cohorts <- unique(myalpha$Cohort)
unique_days    <- unique(myalpha$Day)
levels(unique_days) <- day_levels #set levels of days 

#---------------------------#
####---Alpha Summaries---####
#---------------------------#

# Empty matrix for saving results
range_res <- matrix(NA, nrow = length(myalpha_index), ncol = length(unique_tissues),
                        dimnames = list(myalpha_index, unique_tissues))

## Loop through each tissue and print range of index results
for (tissue in unique_tissues) {
  # Subset the data for the current tissue
  tissue_data <- myalpha %>% filter(Tissue == tissue)
  
  # Loop through each index
  for (index in myalpha_index) {
    # Extract the values for the current index
    index_values <- tissue_data[[index]]
    
    # Calculate the range of values
    index_range <- range(index_values)
    
    # Store the range in the result matrix
    range_res[index, tissue] <- paste(round(index_range, 3), collapse = "-")
  }
}

range_res <- noquote(range_res)
range_res


# Empty lists for saving results
tissue_alpha.MSE <- list()
cohort_alpha.MSE <- list()
group_alpha.MSE  <- list()

## Call summary_stats function
tissue_alpha.MSE <- summary_stats(myalpha, "Tissue", unique_tissues)
cohort_alpha.MSE <- summary_stats(myalpha, "Coh_Tis", unique_CT)
group_alpha.MSE  <- summary_stats(myalpha, "Coh_Tis_Day", unique_CTD)

#----------------------------#
####---Alpha Statistics---####
#----------------------------#

# Shaprio Test
## List of alpha indices of interest
alphaidx <- list(Observed_ASVs = myalpha$observed, 
                 Simpson_Evenness = myalpha$evenness_simpson, 
                 Phylogenetic_Diversity = myalpha$PD)

## For loop Shapiro test, running through indices listed in alphaidx
for(i in seq_along(alphaidx)) {
  index_name <- names(alphaidx)[i]
  test.res <- shapiro.test(alphaidx[[i]])
  
  print(paste0("Index: ", index_name))
  print(test.res) #prints results to screen
}


# Kruskal-Wallis and Dunn's Post-Hoc Tests
## Empty lists for saving results
tissue_alpha_stats <- data.frame() #observed ASVs ~ Cohort for each tissue
cohort_alpha_stats <- data.frame() #observed ASVs ~ Tissue for each cohort

  
## Loop through each cohort
for(tissue in unique_tissues) {
  # Loop through each unique Day
  for(day in unique_days) {
    # Loop through each index
    for(index in myalpha_index) {

    # Subset the data for the current Tissue and Day
    subset_data <- myalpha %>% filter(Tissue == tissue, Day == day)

    # Perform Kruskal-Wallis test
    #result <- kruskal_dunn(index, "Cohort", subset_data)
    result <- wilcox.test(get(index) ~ Cohort, data = subset_data, correct = TRUE, exact = FALSE)

    # Add tissue and day columns to results
    result_df <- data.frame(
      Tissue = tissue,
      Day = day,
      Index = index,
      Statistic = result$statistic,
      P.Value = result$p.value)

    # Bind the result to the results data frame
    cohort_alpha_stats <- rbind(cohort_alpha_stats, result_df)
    }
  }
}


## Loop through each tissue
for(cohort in unique_cohorts) { 
  # Loop through each unique Day
  for (day in unique_days) {
    # Loop through each index
    for(index in myalpha_index) {
      
      # Subset the data for the current Tissue and Day
      subset_data <- myalpha %>% filter(Cohort == cohort, Day == day)
      
      # Perform Kruskal-Wallis test
      result <- kruskal_dunn(index, "Tissue", subset_data)
      
      # Add tissue and day columns to results
      result$Cohort <- cohort
      result$Day <- day
      
      # Bind the result to the results data frame
      tissue_alpha_stats <- rbind(tissue_alpha_stats, result)
    }
  }
}


## Test significance of range differences between hemolymph and hepatopancreas
fligner.test(observed ~ Tissue, data = myalpha.s)

## Clean up tables 
cohort_alpha_stats <- cohort_alpha_stats %>%
  arrange(index)

tissue_alpha_stats <- tissue_alpha_stats %>% 
  mutate(comparison = gsub("Hemolymph", "Hm", comparison),
         comparison = gsub("Hepatopancreas", "Hp", comparison),
         comparison = gsub("Water", "W", comparison)) %>% 
  separate_rows(comparison, Z, p.adj, sep = ",") %>%
  arrange(index)

## Pull only significant results
sig.cohort <- cohort_alpha_stats %>% filter(P.Value <= 0.05) %>% arrange(Index, Tissue, Day)

sig.tissue <- tissue_alpha_stats %>% filter(p.adj <= 0.05) %>% arrange(index, Cohort, Day)

# Save results
write.csv(cohort_alpha_stats, paste0(results_fd, "Cohort_AlphaStats.csv"), row.names = FALSE)
write.csv(tissue_alpha_stats, paste0(results_fd, "Tissue_AlphaStats.csv"), row.names = FALSE)

#-------------------------------#
####---Fixed Model Testing---####
#-------------------------------#

# All tissues included 
f.mod1 <- "Cohort + Tissue + Day" # Additive model
f.mod2 <- "Cohort * Tissue * Day" # Interaction model
alpha_f.model <- c(f.mod1, f.mod2)

# Within tissue models
tf.mod1 <- "Cohort + Day" # Additive model
tf.mod2 <- "Cohort * Day" # Interaction model
alpha_tf.model <- c(tf.mod1, tf.mod2)

# Empty lists for saving results
af.results <- list()    # All tissues included (alpha_f.model)
af.hm.results <- list() # Hemolymph-specific model results (alpha_tf.model)
af.hp.results <- list() # Hepatopancreas-specific model results (alpha_tf.model)

# Call alpha_mod_comp function to compare models
af.results <- alpha_mod_comp(myalpha, myalpha_index, alpha_f.model)
af.hm.results <- alpha_mod_comp(myalpha, myalpha_index, alpha_tf.model, tissue = "Hemolymph")
af.hp.results <- alpha_mod_comp(myalpha, myalpha_index, alpha_tf.model, tissue = "Hepatopancreas")

# Post-hoc analysis for best fitting model
## Hemolymph
obs.mod.posthoc.hm <- fixmodel_posthoc("observed", tf.mod1, "Hemolymph")  
even.mod.posthoc.hm <- fixmodel_posthoc("evenness_simpson", tf.mod1, "Hemolymph")
pd.mod.posthoc.hm <- fixmodel_posthoc("PD", tf.mod1, "Hemolymph")

## Hepatopancreas
obs.mod.posthoc.hp <- fixmodel_posthoc("observed", tf.mod1, "Hepatopancreas")  
even.mod.posthoc.hp <- fixmodel_posthoc("evenness_simpson", tf.mod1, "Hepatopancreas")
pd.mod.posthoc.hp <- fixmodel_posthoc("PD", tf.mod1, "Hepatopancreas")

#---------------------#
####---Plotting----####
#---------------------#

# Publication figures
plot_theme <- theme(axis.text.x = element_text(vjust = 0.5),
                    axis.title = element_text(size = 12, face = "bold"), 
                    panel.grid.major = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(color = "black"))

## Custom labels to shorten tissue names in facet labels
custom_labels <- c("Hemolymph" = "Hm", "Hepatopancreas" = "Hp", "Water" = "W")

# Observed ASVs 
Obs_cohort <-	ggplot(myalpha, aes(y = observed, x = as.factor(Day), color = Cohort)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.25) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +
  scale_color_manual(values = cohort_cols) +
  scale_y_continuous(breaks = c(100, 300, 500)) + 
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "black") +
  plot_theme +
  labs(x = "", y = "Observed ASVs", tag = "A") 

Obs_cohort <- Obs_cohort + facet_grid(Tissue ~ ., labeller = labeller(Tissue = custom_labels)) + 
  theme(strip.text.x = element_text(size = 11, face = "bold"),
        strip.text.y = element_text(size = 11, face = "bold"))

# Simpson Evenness
Simp_cohort <-	ggplot(myalpha, aes(y = evenness_simpson, x = as.factor(Day), color = Cohort)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.25) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +
  scale_color_manual(values = cohort_cols) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "black") +
  plot_theme +
  labs(x = "", y = "Simpson Evenness", tag = "B") 

Simp_cohort <- Simp_cohort + facet_grid(Tissue ~ ., labeller = labeller(Tissue = custom_labels)) + 
  theme(strip.text.x = element_text(size = 11, face = "bold"),
        strip.text.y = element_text(size = 11, face = "bold"))

# Phylogenetic Diversity
PD_cohort <-	ggplot(myalpha, aes(y = PD, x = as.factor(Day), color = Cohort)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.25) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +
  scale_color_manual(values = cohort_cols) +
  scale_y_continuous(breaks = c(10, 30, 50)) + 
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "black") +
  plot_theme +
  labs(x = "Day", y = "Phylogenetic Diversity", tag = "C") 

PD_cohort <- PD_cohort + facet_grid(Tissue ~ ., labeller = labeller(Tissue = custom_labels)) + 
  theme(strip.text.x = element_text(size = 11, face = "bold"), 
        strip.text.y = element_text(size = 11, face = "bold"))

# Save as PDF
pdf(paste0(results_fd, "Figure3.pdf"), width = 8, height = 10)

ggarrange(Obs_cohort, Simp_cohort, PD_cohort, nrow = 3, ncol = 1, common.legend = TRUE, legend = "bottom")

dev.off()
