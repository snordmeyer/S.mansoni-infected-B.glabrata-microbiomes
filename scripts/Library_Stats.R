#-------------------------#
####---Load Packages---####
#-------------------------#

library(readr)     #read_tsv
library(magrittr)  #pipe operators
library(ggplot2)   #plotting
library(ggpubr)    #ggarrange function

#------------------------#
####---Import files---####
#------------------------#

# Folders
QIIME_fd <- "results/0 - QIIME2/"
results_fd <- "results/1 - Library statistics/"

# Sequencing stat file
myqzv.f <- paste0(QIIME_fd, "denoising-stats.qzv")

# Loading data
## Identification of the metadata file within the archive file
myfile <- unzip(myqzv.f, list=TRUE)[1] %>% unlist() %>% grep("metadata.tsv", ., value=T)
## Reading the file from the archive file
libstats.f <- read.delim(unz(myqzv.f, myfile), stringsAsFactors=FALSE, check.names=FALSE)

# libstats.f <- read_tsv(paste0(QIIME_fd, "denoising-metadata.tsv")) #reads info

# Identify numeric columns
num <- which(libstats.f[1,] == "numeric")

# Remove q2:types column information
libstats.f <- libstats.f[-1,]

# Convert values to numeric
libstats.f[ ,num] <- sapply(libstats.f[ ,num], function(x) as.numeric(x))

#-------------------------#
####---Average Reads---####
#-------------------------#
source("scripts/Functions.R")

if(! dir.exists(results_fd)) { dir.create(results_fd, recursive = TRUE) }

# Average number of input reads 
inputreads <- summary(libstats.f$input)
MSE(libstats.f$input)

# Average number of filtered reads 
filteredreads <- summary(libstats.f$filtered)
MSE(libstats.f$filtered)

# Average number of filtered reads by sample type 
## Loop through each unique tissue
unique_tissues <- unique(libstats.f$Tissue)

## Empty list for saving results
tissue_results <- list()

for (tissue in unique_tissues) {
  subset_data <- libstats.f[libstats.f$Tissue == tissue, ]
  result <- MSE(subset_data$`non-chimeric`)  
  tissue_results[[tissue]] <- result
}

## Loop through each unique cohort
unique_cohorts <- unique(libstats.f$Cohort)

## Empty list for saving results
cohort_results <- list()

for (group in unique_cohorts) {
  subset_data <- libstats.f[libstats.f$Cohort == group, ]
  result <- MSE(subset_data$`non-chimeric`)  
  cohort_results[[group]] <- result
}

#---------------------------#
####---Write out table---####
#---------------------------#

mytable <- libstats.f[, c("id", "Tray", "Snail_Number", "Day", "Tissue", "Cohort", "input", "filtered", "denoised", "non-chimeric") ]
colnames(mytable) <- c("Sample ID", "Tray", "Snail_Number", "Day", "Tissue", "Cohort", "Number of input sequences", "Number of sequences after filtering", "Number of sequences after denoising", "Number of non chimeric sequences")

write.table(mytable, paste0(results_fd, "SuppTable1.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

#--------------------#
####---Plotting---####
#--------------------#

# Set plot theme
lib_theme <- theme(axis.text.x = element_text(vjust = 0.5),
                   axis.text = element_text(face = "bold"),
                   axis.title.y = element_text(face = "bold"),
                   axis.title.x = element_blank(),
                   plot.title = element_text(hjust = 0.5),
                   panel.grid.major = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(color = "#101010"))

# Set plot colors
cohort_cols <- c("#2C7800", "#b54e5d")

# Plot reads by type
input_plot <-	ggplot(libstats.f, aes(y = input, x = Tissue, color = Cohort)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +
  scale_color_manual(values = cohort_cols) +
  lib_theme +
  labs(y = "Number of Reads", title = "Input Reads") 

filtered_plot <-	ggplot(libstats.f, aes(y = filtered, x = Tissue, color = Cohort)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +
  scale_color_manual(values = cohort_cols) +
  lib_theme +
  labs(y = "Number of Reads", title = "Filtered Reads") 

denoised_plot <-	ggplot(libstats.f, aes(y = denoised, x = Tissue, color = Cohort)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +
  scale_color_manual(values = cohort_cols) +
  lib_theme +
  labs(y = "Number of Reads", title = "Denoised Reads") 

nonchimeric_plot <-	ggplot(libstats.f, aes(y = `non-chimeric`, x = Tissue, color = Cohort)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +
  scale_color_manual(values = cohort_cols) +
  lib_theme +
  labs(y = "Number of Reads", title = "Non-Chimeric Reads") 

# Print plots as pdf
pdf(paste0(results_fd, "SuppFigure2.pdf"), width = 8, height = 8)

ggarrange(input_plot, filtered_plot, denoised_plot, nonchimeric_plot, 
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
  
dev.off()
