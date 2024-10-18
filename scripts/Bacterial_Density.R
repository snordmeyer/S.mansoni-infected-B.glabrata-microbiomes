#-------------------------#
####---Load Packages---####
#-------------------------#

# File management
library(readxl)
library(writexl)

# Data cleaning
library(tidyverse)
library(tidyr)
library(dplyr)
library(rqdatatable)

# Plotting
library(ggplot2)
library(VennDiagram)

# Statistics
library(car)
library(FSA)
library(rstatix)
library(ggpubr)
library(rcompanion)

#------------------------#
####---Import files---####
#------------------------#

results_fd <- "results/4 - Bacterial density/"
data_fd <- "data/"

masterfile <- read_excel(paste0(data_fd, "Hm_16s_qPCR.xlsx"))

source("scripts/Functions.R")

# Set colors
cohort_cols <- c("#2C7800", "#b54e5d")

# Volumes
elu_vol <- 50
hml_vol <- 40

#----------------------------------#
####---Convert Cq to Quantity---####
#----------------------------------#

# Set slope and y-intercept from qPCR standard curve
slope <- -4.48 
intercept <- 44.719

# Add Quantity to masterfile
masterfile <- masterfile %>%
  mutate(Quantity = estimate_quantity(CT_mean, slope, intercept))

masterfile

#-------------------------------------#
####---Correct and normlize data---####
#-------------------------------------#
# Dilution factor to reflect sample concentration (instead of elution concentration)
dil_fct <- elu_vol/hml_vol

hm <- as.data.frame(masterfile)

# Reorder row
hm <- hm[ order(hm[, "Sample_ID"]), ]

# Normalization
hm$Norm_Quant <- hm[, "Quantity"] * 10 * dil_fct

hm$Cohort <- sapply(hm[, "Tray"], function(x) ifelse(x == "EA" | x == "EB", "Exposed", "Control"))


# Clean data
hm <- hm %>%
  mutate(Day = factor(Day, levels = c("0", "7", "14", "21", "28", "35", "37", "39"))) %>%
  filter(Norm_Quant < 500000)

hm$Cohort <- as.factor(hm$Cohort)
hm$Day <- as.factor(hm$Day)

# Set variables
unique_cohort <- unique(hm$Cohort)
unique_day <- unique(hm$Day)

#-----------------------------#
####---Density summaries---####
#-----------------------------#

if(! dir.exists(results_fd)) { dir.create(results_fd, recursive = TRUE) }

for(group in unique_cohort) {
  
  #subset data
  subset <- hm %>% filter(Cohort == group) %>% summary()
  
  #print results
  print(group)
  print(subset)
  
}


#------------------------#
####---Shaprio Test---####
#------------------------#

hm %>% 
  group_by(Day, Cohort) %>%
  shapiro_test(Norm_Quant) # if p > 0.05, data is normal ; if not, deviates from normal distribution

#---------------------------------------------------#
####---Kruskal-Wallis and Dunn's Post-Hoc Test---####
#---------------------------------------------------#

# Empty lists for saving results
cohort_res <- vector("list")
day_res <- vector("list")

# Loop through each unique Cohort
for(group in unique_cohort) {
  # Subset the data for the current cohort
  subset <- hm %>% filter(Cohort == group)
  
  # Kruskal Wallis test
  kruskal_res <- kruskal.test(Norm_Quant ~ Day, data = subset)
  
  # Dunn's post hoc test
  posthoc <- dunnTest(Norm_Quant ~ Day, data = subset)
  PT <- posthoc$res
  letters <- cldList(comparison = PT$Comparison,
                     p.value = PT$P.adj,
                     threshold = 0.05)
  
  # Save results
  cohort_res[[paste("Cohort", group)]] <- list(Kruskal_Wallis = kruskal_res,
                                               Dunn_PostHoc = posthoc,
                                               PT = PT,
                                               Letters = letters)
  
}

cohort_res

# Loop through each unique Day
for(day in unique_day) {
  # Subset the data for the current day
  subset <- hm %>% filter(Day == day)
  
  # Kruskal Wallis test
  #kruskal_res <- kruskal.test(Norm_Quant ~ Cohort, data = subset)
  result <- wilcox.test(Norm_Quant ~ Cohort, data = subset, correct = TRUE, exact = FALSE)
  
  # Save results
  #day_res[[paste("Day", day)]] <- list(Kruskal_Wallis = kruskal_res)
  day_res[[paste("Day", day)]] <- data.frame(
    Statistic = result$statistic,
    P.Value = result$p.value)
  
}

day_res

# Plotting
## Pulling letters from Dunn's
Control_letters <- cohort_res$`Cohort Control`$Letters[1:2]
Control_letters <- Control_letters %>% mutate(Cohort = "Control") %>% mutate(Group = if_else(Group == "", "0", Group))

Exposed_letters <- cohort_res$`Cohort Exposed`$Letters[1:2]
Exposed_letters <- Exposed_letters %>% mutate(Cohort = "Exposed") %>% mutate(Group = if_else(Group == "", "0", Group))

day_positions <- data.frame(
  Day = levels(hm$Day),  # Get levels of Day
  x_numeric = seq_along(levels(hm$Day))  # Create a sequence for numeric positions
)

dunn_letters <- bind_rows(Control_letters, Exposed_letters)
dunn_letters <- dunn_letters %>% rename(Day = Group)
dunn_letters <- dunn_letters %>%
  left_join(day_positions, by = "Day") %>%  # Join to get the numeric positions
  mutate(x_adjusted = case_when(
    Cohort == "Control" ~ x_numeric - 0.2,  # Shift Control slightly left
    Cohort == "Exposed" ~ x_numeric + 0.2    # Shift Exposed slightly right
  ))

## Plot and save
pdf(paste0(results_fd, "Figure2.pdf"), width = 8, height = 4)

hm_plot <- ggplot(hm, aes(x = Day, y = Norm_Quant, color = Cohort)) +
  geom_boxplot(size = 0.7) +
  scale_y_log10() +
  scale_color_manual(values = cohort_cols) +
  labs(y = "16s copies/μL") +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "black", size = 0.7) + 
  geom_text(data = dunn_letters, aes(x = x_adjusted, y = 1e+06, 
                                     label = Letter, color = Cohort),
            size = 5, fontface = "bold", show.legend = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 15),
        legend.key.size = unit(0.75, 'cm'),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "bottom")

hm_plot

dev.off()

#--------------------------#
####---Average copies---####
#--------------------------#

# Empty list for saving results
copies <- vector("list")

# Loop through each cohort
for(group in unique_cohort) {
	# Subset the data for the current cohort
	subset <- hm %>% filter(Cohort == group) %>% na.omit()
	
	# Use MSE function for mean and standard error
	average <- MSE(subset$Norm_Quant)
	
	# Save results
	copies[[paste("Cohort", group)]] <- list(Cohort = group, Average = average)
}

copies
