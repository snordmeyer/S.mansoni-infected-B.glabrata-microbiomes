#-------------------------#
####---Load Packages---####
#-------------------------#

library(dplyr)
library(survival)
library(survminer)
library(readxl)
library(gridExtra)

#------------------------#
####---Import files---####
#------------------------#

data_fd <- "data/"

results_fd <- "results/2 - Survival/"
if(! dir.exists(results_fd)) { dir.create(results_fd, recursive = TRUE) }

surv_data <- read_xlsx(paste0(data_fd, "Aim1_SurvivalData.xlsx")) # only snails not assigned to be sampled

#-------------------------#
####---Survival plot---####
#-------------------------#

fit <- survfit(Surv(Time, Status) ~ Tray, data = surv_data)

fit.summ <- summary(fit)$table

surv_theme <- theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14))

# Print plot as pdf
pdf(paste0(results_fd, "SuppFigure1.pdf"), width = 8, height = 10)

surv_plot <- ggsurvplot(fit, conf.int = FALSE, pval = TRUE, risk.table = FALSE,
              legend.labs=c("CA", "CB", "EA", "EB"), legend.title = "Tray",
              title = "Kaplan-Meier Survival Curve",
              xlab = "Time (Days)", ylab = "Survival",
              risk.table.height = 0.25, 
              ggtheme = surv_theme,
              palette = c("#2C7800","#5FA052","#b54e5d","#b59096"))

grid.arrange(surv_plot$plot, ncol = 1)

dev.off()

#-------------------------#
####---Survival rate---####
#-------------------------#

Trays <- c("CA", "CB", "EA", "EB")
Survived <- c(11, 15, 3, 2)
Not_Sampled_Total <- c(22, 22, 24, 24)

Survival_Calc <- data.frame(Trays, Survived, Not_Sampled_Total)

Survival_Calc <- Survival_Calc %>% 
  mutate(Survival_Rate = Survived/Not_Sampled_Total * 100, 
         Mortality_Rate = (Not_Sampled_Total-Survived)/Not_Sampled_Total * 100)

#-------------------------#
####---Log-rank test---####
#-------------------------#

surv_diff <- survdiff(Surv(Time, Status) ~ Cohort, 
                      data = surv_data)
surv_diff

# Control only 
surv_data_C <- surv_data %>% filter(Cohort == "Control")
surv_diff_C <- survdiff(Surv(Time, Status) ~ Tray, 
                        data = surv_data_C)
surv_diff_C

# Exposed only 
surv_data_E <- surv_data %>% filter(Cohort == "Exposed")
surv_diff_E <- survdiff(Surv(Time, Status) ~ Tray, 
                        data = surv_data_E)
surv_diff_E

#----------------------------------#
####---Pairwise post-hoc test---####
#----------------------------------#

pairwise_diff <- pairwise_survdiff(Surv (Time, Status) ~ Cohort, 
                                   data = surv_data, p.adjust.method = "bonferroni") 
pairwise_diff

# Control only 
pairwise_diff_C <- pairwise_survdiff(Surv (Time, Status) ~ Tray, 
                                     data = surv_data_C, p.adjust.method = "bonferroni") 
pairwise_diff_C

# Exposed only 
pairwise_diff_E <- pairwise_survdiff(Surv (Time, Status) ~ Tray, 
                                     data = surv_data_E, p.adjust.method = "bonferroni") 
pairwise_diff_E

