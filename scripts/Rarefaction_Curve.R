#-------------------------#
####---Load Packages---####
#-------------------------#

library(ranacapa)
library(ggplot2)

#-----------------------------#
####---Rarefaction Curve---####
#-----------------------------#

# Folders
data_fd <- "results/3 - Create Phyloseq/"

results_fd <- "results/1 - Library statistics/"
if(! dir.exists(results_fd)) { dir.create(results_fd, recursive = TRUE) }

# Load mybiom2 object
mybiom2 <- readRDS(paste0(data_fd, "mybiom2.rds"))

# Set plot colors
cohort_cols <- c("#2C7800", "#b54e5d")

# Print plots as pdf
pdf(paste0(results_fd, "Figure3.pdf"), width = 8, height = 8)

rc <- ggrare(mybiom2, step = 1000, color = "Cohort", se = FALSE, plot=FALSE) +
  geom_line(size = 0.75) +
  scale_color_manual(values = cohort_cols) +
  theme(legend.position="bottom",
        legend.key = element_blank(), legend.key.width = unit(4,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        strip.text.x = element_text(size = 11, face = "bold"),
        strip.text.y = element_text(size = 11, face = "bold")) +
  facet_grid(Tissue ~ Day)

rc

dev.off()
