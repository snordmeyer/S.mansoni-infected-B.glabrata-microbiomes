#-------------------------#
####---Load Packages---####
#-------------------------#

library(tidyr)
library(microViz)
library(phyloseq)
library(phyloseq.extended)
library(scales)
library(ape)
library(ggplot2)
library(speedyseq)

#---------------------#
####---Load Data---####
#---------------------#

data_fd <- "results/3 - Create Phyloseq/"

results_fd <- "results/7 - Taxonomy/"

mybiom2 <- readRDS(paste0(data_fd, "mybiom2.rds"))

# Set levels
day_levels <- c("0", "7", "14", "21", "28", "35", "37", "39")

cohort_cols <- c("#2c7800", "#b54e5d")

#-----------------------#
####---Aggregation---####
#-----------------------#

# Aggregate at rank = "phylum"
## Rearrange tax_table(mybiom2) layout - otherwise will have errors
tax_table(mybiom2) <- tax_table(mybiom2)[,c(2:9,1)]

## Create new mybiom for rank filtering
mybiom3 <- mybiom2 %>%
  tax_fix() %>%
  tax_agg(rank = "Phylum")


#-----------------------------------------------#
####---Composition relative abundance plot---####
#-----------------------------------------------#

# Create merged data
mybiom_gr <- mutate_sample_data(mybiom3, MergeCrit = (concated_column = paste(Cohort, Tissue, Day)))

## Merge samples using MergeCrit - save as new biom object
mymerge <- merge_samples(mybiom_gr, "MergeCrit")

# Create new sample_data object (empty)
col_list <- c("ID", "Cohort", "Tissue", "Day")
mergedSD <- data.frame(matrix(nrow = 48, ncol = length(col_list)))
colnames(mergedSD) <- col_list

## Add in sample_data to mergedSD
mergedSD <- mergedSD %>%
  dplyr::mutate(ID = paste(rownames(sample_data(mymerge)))) %>%
  separate(ID, c("Cohort", "Tissue", "Day"), sep = " ", remove = FALSE)

# Combine mymerge with mergedSD data
OTU_tbl <- otu_table(mymerge, taxa_are_rows = FALSE)
tax_tbl <- tax_table(mymerge)
samps <- sample_data(mergedSD)
sample_names(samps) <- mergedSD$ID

mergedphylo <- phyloseq(OTU_tbl, tax_tbl, samps)

## Create new column - merge of Tissue & Day only (will try to facet by Cohort)
mergedphylo <- mutate_sample_data(mergedphylo, TissDay = (concated_column = paste(Tissue, Day)))

## Change Day from character to numeric & set levels
sample_data(mergedphylo)$Day <- as.numeric(sample_data(mergedphylo)$Day)
levels(sample_data(mergedphylo)$Day) <- day_levels

## Set tissue levels (reverse order for correct plotting)
sample_data(mergedphylo)$Tissue <- factor(sample_data(mergedphylo)$Tissue,
                                          levels = c("Water", "Hepatopancreas", "Hemolymph"))

# Pull specific phyla for histograms
target_phyla <- c("Fusobacteria", "Actinobacteria", "Acidobacteria", "Nitrospirae", "Cyanobacteria", "Firmicutes")

## tax_table to identify taxa in target_phyla
tax <- tax_table(mergedphylo)
matching_taxa <- rownames(tax)[tax[, "Phylum"] %in% target_phyla]

## Prune new phyloseq object
phy_df <- prune_taxa(matching_taxa, mergedphylo)

df_melt <- psmelt(phy_df)

df_melt$Day <- as.factor(df_melt$Day)
df_melt$Tissue <- factor(df_melt$Tissue, levels = c("Hemolymph", "Hepatopancreas", "Water"))


#--------------------#
####---Plotting---####
#--------------------#

if(! dir.exists(results_fd)) { dir.create(results_fd, recursive = TRUE) }
pdf(paste0(results_fd, "Figure5.pdf"), width = 8.5, height = 10)

mergedphylo %>%
  comp_barplot(tax_level = "Phylum", n_taxa = 19, merge_other = TRUE, x = "Tissue", bar_outline_color = NA) +
  facet_grid(cols = vars(Cohort), rows = vars(Day),
             scales = "free", space = "free_y") +
  coord_flip() +
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        strip.background = element_rect(color = "black", fill = "grey"))+
  labs(x = " ")

dev.off()

# Specific phyla barplots

pdf(paste0(results_fd, "SuppFigure6.pdf"), width = 6, height = 8)

ggplot(df_melt, aes(x = Day, y = Abundance, fill = Cohort)) +
	geom_col(position = "dodge") +
	facet_grid(Phylum ~ Tissue, scales = "free_y") +
	theme_minimal() +
	theme(
		strip.text = element_text(face = "bold"),
		panel.border = element_blank(),
		legend.position = "bottom",
		strip.background = element_blank()) +
	labs(y = "Abundance", x = "Day") +
	scale_fill_manual(values = cohort_cols)

dev.off()
