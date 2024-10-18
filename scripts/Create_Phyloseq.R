#-------------------------#
####---Load Packages---####
#-------------------------#

library(qiime2R)
library(magrittr)
library(tidyr)
library(dplyr) 
library(stringr)
library(phyloseq)
library(speedyseq)

#------------------------#
####---Import files---####
#------------------------#

QIIME_fd <- "results/0 - QIIME2/"
results_fd <- "results/3 - Create Phyloseq/"

# Load metadata
md <- read.delim("data/A1_metadata.txt", sep = '\t', header = TRUE)
row.names(md) <- md$sample.id

# Load QIIME files
aligned.f <- qiime2R::read_qza(paste0(QIIME_fd, "aligned-rep-seqs.qza"))
table.f <- qiime2R::read_qza(paste0(QIIME_fd, "table.qza"))
taxa.f <- qiime2R::read_qza(paste0(QIIME_fd, "rep-seqs_taxa.qza"))
unass.f <- qiime2R::read_qza(paste0(QIIME_fd, "rep-seqs_unassigned.qza"))
tree.f <- qiime2R::read_qza(paste0(QIIME_fd, "rooted-tree.qza"))

blast <- read.delim(paste0(QIIME_fd, "rep-seqs-unassigned/unassigned.blastn.tsv"), header = FALSE)

#-----------------------------------------#
####---Create variables & clean data---####
#-----------------------------------------#

asv <- as.data.frame(table.f$data)
asv.tb <- (otu_table(asv, taxa_are_rows = TRUE))

taxa <- as.data.frame(taxa.f$data)

tree <- tree.f$data

## Splitting information
cln.nm <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa <- separate(taxa, 2, cln.nm, sep = ";")

## Cleaning names
taxa[, (1:length(cln.nm))+1] <- apply(taxa[,(1:length(cln.nm))+1], 2, function(x) gsub(".*__", "", x))

## Polishing data
rownames(taxa) <- taxa[,1]
taxa[,-c(1,ncol(taxa))]

## Assign "unassigned" to all levels of unassigned ASVs
taxa[ taxa[,2] == "Unassigned", ] <- "Unassigned"

## Convert table
taxa.tb <- phyloseq::tax_table(as.matrix(taxa))
taxa_names(taxa.tb) <- rownames(taxa)

#------------------------------------#
####---Building phyloseq object---####
#------------------------------------#

mybiom <- phyloseq(asv.tb, taxa.tb, tree, sample_data(md))

## Cleaning data from contaminants
nb.asv <- nrow(phyloseq::tax_table(mybiom))
nb.unass <- (taxa[,1] == "Unassigned") %>% sum()

# Mitochondira & chloroplast contaminants
my.mt <- apply(phyloseq::tax_table(mybiom), 1, function(x) any(grepl("mitochondria", x, ignore.case = TRUE)))
my.cp <- apply(phyloseq::tax_table(mybiom), 1, function(x) any(grepl("chloroplast", x, ignore.case = TRUE)))

# Eukaryote contaminants 
my.euk <- blast[ blast[,14] == "Eukaryota", 1]
my.euk <- rownames(phyloseq::tax_table(mybiom)) %in% my.euk 

# Exporting contaminant list
write.table(rownames(phyloseq::tax_table(mybiom))[(my.mt | my.cp | my.euk)], "data/contaminants", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Remove contaminants
mybiom2 <- subset_taxa(mybiom, ! (my.mt | my.cp))

# Final numbers
cat("\nNumber of ASVs removed:", "\n\t- mitochondria and chloroplasts:", sum(my.mt | my.cp), "/", signif(sum(my.mt | my.cp) / nb.asv * 100, 2), 
    "% of the total ASVs:", "\n\t- eukaryotes:", sum(my.euk), "/", signif(sum(my.euk) / nb.unass * 100, 2), "% of the unassigned ASVs", "/", 
    signif(sum(my.euk) / nb.asv * 100, 2), "% of the total ASVs", "\n\t- total:", sum(my.mt | my.cp | my.euk), "/", signif(sum(my.mt | my.cp | my.euk) / nb.asv *100, 2),
    "% of the total ASVs", "/n/n")

## Number of ASVs after filtering
filtered.asv <- nb.asv - sum(my.mt | my.cp | my.euk)

#----------------------------------#
####---Adjust phyloseq object---####
#----------------------------------#

# Create new columns for cohort grouping
new_cols <- md %>%
  mutate(Day = replace(Day, Day == 1, 0),
         sample.id = str_replace(sample.id, "^(.*?\\-\\d\\-)1(\\-.*)$", "\\10\\2"),
         Coh_Tis_Day = (concated_column = paste(Cohort, Tissue, Day)),
         Coh_Tis = (concated_column = paste(Cohort, Tissue)),
         Coh_Day = (concated_column = paste(Cohort, Day)),
         Tis_Day = (concated_column = paste(Tissue, Day))) %>%
  unite(SnailID, c(Tray, Snail_Number, Day), sep = "_", remove = FALSE) 

## Add new_cols to mybiom2
sample_data(mybiom2) <- new_cols

# Set levels for mybiom2
day_levels <- c("0", "7", "14", "21", "28", "35", "37", "39")
tissue_levels <- c("Hemolymph", "Hepatopancreas", "Water")

# Set "Day" as factor
sample_data(mybiom2)$Day <- as.factor(sample_data(mybiom2)$Day)

## Set levels
levels(sample_data(mybiom2)$Day) <- day_levels

# Save mybiom2 object
if(! dir.exists(results_fd)) { dir.create(results_fd, recursive = TRUE) }
saveRDS(mybiom2, file = paste0(results_fd, "mybiom2.rds"))
saveRDS(tree, file = paste0(results_fd, "tree.rds"))
saveRDS(asv, file = paste0(results_fd, "asv.rds"))
saveRDS(new_cols, file = paste0(results_fd, "new_cols.rds"))
