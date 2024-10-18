#---------------------------------#
####---- Bacterial Density ----####
#---------------------------------#

####---CT to Quantity---####
#--------------------------#
estimate_quantity <- function(CT_mean, slope, intercept) {
  Quantity <- 10^((CT_mean - intercept) / slope)
  return(Quantity)
}

#-------------------------------#
####---- Alpha Diversity ----####
#-------------------------------#

####---Calculating Mean and Standard Error---####
#-----------------------------------------------#

MSE <- function(data) {
  mean_value <- mean(data)
  se_value <- sd(data)/sqrt(length(data))
  return(c(mean = mean_value, se = se_value))
}

####----Calling MSE function over each index---####
#-------------------------------------------------#

summary_stats <- function(myalpha, grouping_var, unique_values) {
  # Empty list for saving results
  results <- list()
  
  for (variable in unique_values) {
    # Subset data by current grouping variable
    subset_data <- myalpha %>% filter(.data[[grouping_var]] == variable)
    
    # Use MSE function on each alpha index
    obs.result <- MSE(subset_data$observed)
    even.result <- MSE(subset_data$evenness_simpson)
    PD.result <- MSE(subset_data$PD)
    
    # Save results in list
    results[[variable]] <- list(Level = variable,
                                Obs.ASV = obs.result, 
                                Simpson = even.result,
                                PD = PD.result)
  }
  
  return(results)
}

####---Kruskal and Dunn's tests---####
#------------------------------------#

# Declare 'index' (observed, evenness_simpson, PD) and 'variable' for comparison (~ Cohort or ~ Tissue)

kruskal_dunn <- function(index, variable, data_set) {

  # Table for saving results
  results_df <<- data.frame(index = character(),
                            x2 = numeric(),
                            df = numeric(),
                            p.value = numeric(),
                            Z = character(),
                            p.adj = character(),
                            comparison = character(),
                            stringsAsFactors = FALSE)

  # Kruskal & Dunn's tests
  ## Set formula
  formula_str <- paste0(index, " ~ as.factor(", variable, ")", collapse = " ")
  formula <- as.formula(formula_str)

  ## Run tests
  kruskal.res <- kruskal.test(formula, data = subset_data)
  dunn.res <- dunnTest(formula, data = subset_data, method = "bonferroni")

  # Save results
  results_df[nrow(results_df) + 1, ] <<- c(index = index,
                                           x2 = round(kruskal.res$statistic, 5),
                                           df = kruskal.res$parameter,
                                           p.value = round(kruskal.res$p.value, 5),
                                           Z = paste(round(dunn.res$res$Z, 5), collapse = ","),
                                           p.adj = paste(round(dunn.res$res$P.adj, 5), collapse = ","),
                                           comparison = paste(dunn.res$res$Comparison, collapse = ","),
                                           stringsAsFactors = FALSE)

  return(results_df)

}

####---Alpha model comparison---####
#----------------------------------#

alpha_mod_comp <- function(data, index, alpha_models, tissue = NULL) {
  # Lists for saving results
  results <- vector("list", length(index))
  
  # Loop through each index variable
  for (idx in index) {
    # Set data subset based on tissue
    subset_data <- if (is.null(tissue)) {
      data  # Include all tissues
    } else {
      data %>% filter(Tissue == tissue)  # Filter by specified tissue
    }
    
    # Make sure Day is set to factor
    subset_data$Day <- as.factor(subset_data$Day)
    
    # Initialize results for the current index
    results[[idx]] <- vector("list", length(alpha_models))
    
    # Test each model
    for (model in alpha_models) {
      results[[idx]][[model]] <- list()
      
      # Run model
      mod <- lm(as.formula(paste(idx, "~", model)), data = subset_data, REML = FALSE) # REML to FALSE because doing model selection
      
      # Model results
      res_bic <- BIC(mod)
      res_aic <- AIC(mod)
      ano <- anova(mod)
      ansum <- summary(ano)
      
      # Save results
      results[[idx]][[model]] <- list(
        index = idx,
        model = mod,
        BIC_score = res_bic,
        AIC_score = res_aic,
        anova = ano,
        summary = ansum
      )
    }
  }
  
  return(results)
}


####---Model Post-Hoc analysis---####
#-----------------------------------#

fixmodel_posthoc <- function(index, model, tissue) {
  # Set data subset
  subset_data <- myalpha %>% filter(Tissue == tissue)
  
  # Make sure Day is set to factor
  subset_data$Day <- as.factor(subset_data$Day)
  
  # Run model
  mod <- lm(paste(index, "~", model, collapse = " "), data = subset_data) #REML to FALSE because doing model selection 
  
  # Post-hoc analysis
  ## Cohort 
  marg.coh <- lsmeans(mod, ~ Cohort)
  CLD.coh <- multcomp::cld(marg.coh, alpha = 0.05, Letters = letters)
  pair.coh <- pairs(marg.coh)
  
  ## Day 
  marg.day <- lsmeans(mod, ~ Day)
  CLD.day <- multcomp::cld(marg.day, alpha = 0.05, Letters = letters)
  pair.day <- pairs(marg.day)
  
  ## Nested
  marginal <- lsmeans(mod, ~ Cohort | Day, options = list(estName = "pred.breaks"))
  CLD.nest <- multcomp::cld(marginal, alpha = 0.05, Letters = letters)
  pair.nest <- pairs(marginal)
  
  results <- list(Cohort_Letters = CLD.coh, Cohort_Padj = pair.coh,
                     Day_Letters = CLD.day, Day_Padj = pair.day,
                     Nested_Letters = CLD.nest, Nested_Padj = pair.nest)
  
  return(results)
}

#------------------------------#
####---- Beta Diversity ----####
#------------------------------#

####---Plotting P and F values in matrix---####
#---------------------------------------------#

plotPFvalues <- function(pvalue.matrix, alg.order=NULL, 
                         show.pvalue=TRUE, font.size=NULL, 
                         border.col="white", pval.col="black", 
                         scientific = TRUE, mybreaks, plot.title=NULL) {
  
  # Dependencies
  if (!requireNamespace("ggplot2", quietly=TRUE)) {
    stop("This function requires the ggplot2 package. Please install it.", call.=FALSE)
  }
  suppressMessages(library("reshape2"))
  
  # Convert the matrix into a data frame and order the algorithms 
  df <- melt(pvalue.matrix)
  colnames(df) <- c("X", "Y", "value")
  df$X <- factor(df$X, levels=day_levels)
  df$Y <- factor(df$Y, levels=day_levels)
  
  ## Create subset of palette 
  mypalette <- RColorBrewer::brewer.pal(4, "Blues")
  
  gplot <- ggplot(df, aes(x = X, y = Y, fill = value)) +
    geom_text(aes(label = scales::number(value)), color = "transparent") +
    geom_tile(data = ~ subset(.x, as.numeric(X) > as.numeric(Y)), aes(fill = value)) + 
    scale_fill_gradientn(colors = mypalette, name = "F value") +
    # Declare new fill scale for the second layer
    new_scale_fill() +
    guides(fill = guide_colorbar(title = "P value")) +
    geom_tile(data = ~ subset(.x, as.numeric(X) < as.numeric(Y)), aes(fill = value)) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "white",
                         limits = c(NA, 0.05), midpoint = 0.05, guide = "colorbar", na.value = "white",
                         breaks = mybreaks, 
                         label = format(mybreaks, digit = 2, scientific = FALSE)) + 
    coord_fixed() +
    theme(panel.background = element_rect(fill = "lightgrey"),
          panel.grid = element_blank(),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 16)) +
    ggplot2::labs(x="F-value" , y="P-value", title=plot.title)
  
  if (show.pvalue) {
    ## geom_text size are not font size
    ## source :https://community.rstudio.com/t/why-does-ggplot-size-parameter-not-behave-consistently/21619/4
    if (is.null(font.size)) font.size <- GeomLabel$default_aes$size / .pt
    p.value.f <- df$value
    if (scientific) {
      p.value.f[ ! is.na(p.value.f) ] <- format(p.value.f[ ! is.na(p.value.f) ], digits=2, 
                                                scientific = FALSE, na.encode = FALSE)
    } else {
      p.value.f <- round(p.value.f, 2)
    }
    gplot <- gplot + ggplot2::geom_text(ggplot2::aes(label = p.value.f),
                                        size=font.size, col=pval.col)
  }
  return(gplot)
}


####---Plot betadisper results in matrix---####
#---------------------------------------------#

plotBDperm <- function(pvalue.matrix, alg.order=NULL, 
                       show.pvalue=TRUE, font.size=NULL, 
                       border.col="white", pval.col="black", 
                       scientific = FALSE, mybreaks, groups) {
  
  # Dependencies
  if (!requireNamespace("ggplot2", quietly=TRUE)) {
    stop("This function requires the ggplot2 package. Please install it.", call.=FALSE)
  }
  suppressMessages(library("reshape2"))
  
  # Convert the matrix into a data frame and order the algorithms 
  df <- melt(pvalue.matrix) # pvalue.matrix
  colnames(df) <- c("X", "Y", "p.value")
  df$X <- factor(df$X, levels=groups)
  df$Y <- factor(df$Y, levels=groups)
  
  # gplot <- ggplot2::ggplot(df, ggplot2::aes(x=X, y=Y, fill=p.value)) + 
  #   ggplot2::geom_tile(col=border.col) +
  #   ggplot2::scale_fill_continuous("p-value") + 
  #   ggplot2::labs(x="Algorithm" , y="Algorithm")
  
  
  gplot <- ggplot(df, aes(x = X, y = Y, fill = p.value)) +
    geom_text(aes(label = scales::number(p.value)), color = "transparent") +
    geom_tile(data = ~ subset(.x, as.numeric(X) > as.numeric(Y)), aes(fill = p.value)) + 
    scale_fill_gradient2(low = "red", mid = "yellow", high = "white",
                         limits = c(NA, 0.05), midpoint = 0.05, guide = "colorbar", na.value = "white",
                         breaks = mybreaks, 
                         label = format(mybreaks, digit = 2, scientific = FALSE)) +
    # Declare new fill scale for the second layer
    new_scale_fill() +
    guides(fill = guide_colorbar(title = "P value")) +
    geom_tile(data = ~ subset(.x, as.numeric(X) < as.numeric(Y)), aes(fill = p.value)) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "white",
                         limits = c(NA, 0.05), midpoint = 0.05, guide = "colorbar", na.value = "white",
                         breaks = mybreaks, 
                         label = format(mybreaks, digit = 2, scientific = FALSE)) + 
    coord_fixed() +
    theme(panel.background = element_rect(fill = "lightgrey"),
          panel.grid = element_blank(),
          axis.text.x = element_text(face = "bold", angle = 90, size = 15),
          axis.text.y = element_text(face = "bold", size = 15),
          axis.title = element_text(size = 16),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          legend.position = "none")
  
  if (show.pvalue) {
    ## geom_text size are not font size
    ## source :https://community.rstudio.com/t/why-does-ggplot-size-parameter-not-behave-consistently/21619/4
    if (is.null(font.size)) font.size <- GeomLabel$default_aes$size / .pt
    p.value.f <- df$p.value
    if (scientific) {
      p.value.f[ ! is.na(p.value.f) ] <- format(p.value.f[ ! is.na(p.value.f) ], digits=2, scientific = FALSE, na.encode = FALSE)
    } else {
      p.value.f <- round(p.value.f, 2)
    }
    gplot <- gplot + ggplot2::geom_text(ggplot2::aes(label = p.value.f),
                                        size=font.size, col=pval.col)
  }
  return(gplot)
}


####---Beta-diversity model analysis---####
#-----------------------------------------#

beta_model <- function(tissue, distance_method, variable, mybiom_data) {        
  
  ## Create distance group
  dist_gr <- mybiom_data
  
  ### Make sure Day is set as factor
  sample_data(dist_gr)[, "Day"] <- as.factor(get_variable(sample_data(dist_gr))[, "Day"]) 
  
  ## Create distance matrix
  dist_mat <- distance(dist_gr, method = distance_method)
  
  ### Basic model testing to identity variables of interest (what may be making an impact?)
  set.seed(my_seed)
  Model <- adonis2(dist_mat ~ Day + Cohort, data = get_variable(sample_data(dist_gr)), permutations = 1000)
  
  # PairwiseAdonis on PERMANOVA results
  set.seed(my_seed)
  Model_post <- pairwise.adonis2(dist_mat ~ Day + Cohort, data = get_variable(sample_data(dist_gr)))
  
  Model_post <- Model_post[c("parent_call", "0_vs_7", "0_vs_14", "0_vs_21", "0_vs_28", "0_vs_35", "0_vs_37", "0_vs_39", 
                           "14_vs_7", "21_vs_7", "28_vs_7", "35_vs_7", "37_vs_7", "39_vs_7", 
                           "14_vs_21","14_vs_28", "14_vs_35", "14_vs_37", "14_vs_39", 
                           "21_vs_28", "21_vs_35", "21_vs_37", "21_vs_39", 
                           "28_vs_35", "28_vs_37", "28_vs_39",
                           "35_vs_37", "35_vs_39", "37_vs_39")]
  
  # Model Comparisons
  Model_AIC <- AICc_permanova2(Model)
  
  ## P-values matrix
  ### Create empty matrices 
  mymat <- matrix(NA, ncol=(length(day_levels)) , nrow=(length(day_levels)))
  colnames(mymat) <- rownames(mymat) <- day_levels
  
  ### Fill in matrix
  if (variable == "Cohort") {
    mymat[lower.tri(mymat, diag = FALSE)] <- as.matrix(as.numeric(unlist(lapply(Model_post[2:29], function(x) x['Pr(>F)'][[1]][2]))))
    mymat <- t(mymat)
    mymat[lower.tri(mymat, diag = FALSE)] <- as.matrix(as.numeric(unlist(lapply(Model_post[2:29], function(x) x['F'][[1]][2]))))
  } else if (variable == "Day") {
    mymat[lower.tri(mymat, diag = FALSE)] <- as.matrix(as.numeric(unlist(lapply(Model_post[2:29], function(x) x['Pr(>F)'][[1]][1]))))
    mymat <- t(mymat)
    mymat[lower.tri(mymat, diag = FALSE)] <- as.matrix(as.numeric(unlist(lapply(Model_post[2:29], function(x) x['F'][[1]][1]))))
  }
  
  # Set breaks
  mybreaks <- c(10^(seq(log10(min(mymat, na.rm=T)), log10(0.05), length.out=3)), 1)
  mybreaks <- trunc(mybreaks*10^3)/10^3 # Cut off decimals to 3 spaces
  
  ### Plot heatmaps
  p1 <- plotPFvalues(mymat, mybreaks = mybreaks, 
               font.size = 4, plot.title = paste0(tissue, ": ", distance_method, " ~ ", variable))

    # Results
  results <- list(
    Tissue = tissue,
    Method = distance_method,
    Model_AIC = Model_AIC,
    Model = Model,
    Model_Post = Model_post,
    Plot = p1
  )
  
  return(results)
  
}




####---Betadisper analysis---####
#-------------------------------#

beta_distance <- function(tissue, distance_method, mybiom_data) {        
  
  ## Create distance group
  dist_gr <- mybiom_data
  sample_data(dist_gr)$Day <- factor(sample_data(dist_gr)$Day)
  
  ## Create distance matrix
  dist_mat <- distance(dist_gr, method = distance_method)
  
  ### Model
  set.seed(my_seed)
  Mod <- adonis2(dist_mat ~ Day + Cohort, data = get_variable(sample_data(dist_gr)), permutations = 1000)
  
  #PairwiseAdonis 
  set.seed(my_seed)
  Mod_post <- pairwise.adonis2(dist_mat ~ Coh_Day, data = get_variable(sample_data(dist_gr)))
  Mod_post <- Mod_post[c("parent_call", 
                         "Control 0_vs_Control 7", "Control 0_vs_Control 14", "Control 0_vs_Control 21", "Control 0_vs_Control 28", "Control 0_vs_Control 35",  
                         "Control 0_vs_Control 37", "Control 0_vs_Control 39", "Control 0_vs_Exposed 0", "Control 0_vs_Exposed 7", "Control 0_vs_Exposed 14","Control 0_vs_Exposed 21", "Control 0_vs_Exposed 28", "Control 0_vs_Exposed 35",  
                         "Control 0_vs_Exposed 37", "Control 0_vs_Exposed 39", "Control 14_vs_Control 7", "Control 21_vs_Control 7", "Control 28_vs_Control 7", "Control 35_vs_Control 7", "Control 37_vs_Control 7", "Control 39_vs_Control 7",
                         "Control 7_vs_Exposed 0", "Control 7_vs_Exposed 7", "Control 7_vs_Exposed 14", "Control 7_vs_Exposed 21", "Control 7_vs_Exposed 28", "Control 7_vs_Exposed 35", "Control 7_vs_Exposed 37", "Control 7_vs_Exposed 39", 
                         "Control 14_vs_Control 21", "Control 14_vs_Control 28", "Control 14_vs_Control 35", "Control 14_vs_Control 37", "Control 14_vs_Control 39", "Control 14_vs_Exposed 0", "Control 14_vs_Exposed 7",  "Control 14_vs_Exposed 14", 
                         "Control 14_vs_Exposed 21", "Control 14_vs_Exposed 28", "Control 14_vs_Exposed 35", "Control 14_vs_Exposed 37", "Control 14_vs_Exposed 39", "Control 21_vs_Control 28", "Control 21_vs_Control 35", "Control 21_vs_Control 37", 
                         "Control 21_vs_Control 39", "Control 21_vs_Exposed 0", "Control 21_vs_Exposed 7", "Control 21_vs_Exposed 14", "Control 21_vs_Exposed 21", "Control 21_vs_Exposed 28", "Control 21_vs_Exposed 35", "Control 21_vs_Exposed 37", 
                         "Control 21_vs_Exposed 39", "Control 28_vs_Control 35", "Control 28_vs_Control 37", "Control 28_vs_Control 39", "Control 28_vs_Exposed 0",  "Control 28_vs_Exposed 7", "Control 28_vs_Exposed 14", "Control 28_vs_Exposed 21",
                         "Control 28_vs_Exposed 28", "Control 28_vs_Exposed 35", "Control 28_vs_Exposed 37", "Control 28_vs_Exposed 39", "Control 35_vs_Control 37", "Control 35_vs_Control 39", "Control 35_vs_Exposed 0", "Control 35_vs_Exposed 7", 
                         "Control 35_vs_Exposed 14", "Control 35_vs_Exposed 21", "Control 35_vs_Exposed 28", "Control 35_vs_Exposed 35", "Control 35_vs_Exposed 37", "Control 35_vs_Exposed 39", "Control 37_vs_Control 39", "Control 37_vs_Exposed 0", 
                         "Control 37_vs_Exposed 7", "Control 37_vs_Exposed 14", "Control 37_vs_Exposed 21", "Control 37_vs_Exposed 28", "Control 37_vs_Exposed 35", "Control 37_vs_Exposed 37", "Control 37_vs_Exposed 39", "Control 39_vs_Exposed 0", 
                         "Control 39_vs_Exposed 7", "Control 39_vs_Exposed 14", "Control 39_vs_Exposed 21", "Control 39_vs_Exposed 28", "Control 39_vs_Exposed 35", "Control 39_vs_Exposed 37", "Control 39_vs_Exposed 39", "Exposed 0_vs_Exposed 7",  
                         "Exposed 0_vs_Exposed 14", "Exposed 0_vs_Exposed 21", "Exposed 0_vs_Exposed 28", "Exposed 0_vs_Exposed 35", "Exposed 0_vs_Exposed 37", "Exposed 0_vs_Exposed 39", "Exposed 14_vs_Exposed 7", "Exposed 21_vs_Exposed 7",  
                         "Exposed 28_vs_Exposed 7", "Exposed 35_vs_Exposed 7", "Exposed 37_vs_Exposed 7", "Exposed 39_vs_Exposed 7","Exposed 14_vs_Exposed 21", "Exposed 14_vs_Exposed 28", "Exposed 14_vs_Exposed 35", "Exposed 14_vs_Exposed 37", 
                         "Exposed 14_vs_Exposed 39", "Exposed 21_vs_Exposed 28", "Exposed 21_vs_Exposed 35", "Exposed 21_vs_Exposed 37", "Exposed 21_vs_Exposed 39", "Exposed 28_vs_Exposed 35", "Exposed 28_vs_Exposed 37", "Exposed 28_vs_Exposed 39",
                         "Exposed 35_vs_Exposed 37", "Exposed 35_vs_Exposed 39", "Exposed 37_vs_Exposed 39")]
  
  # beta-disper test on each of the variables
  bd <- betadisper(dist_mat, type = "centroid", group = get_variable(sample_data(dist_gr))[, "Coh_Day"])
  tukey <- TukeyHSD(bd)  #post-hoc
  ano <- anova(bd)       #anova
  permu <- permutest(bd) #permutest
  
  # Plotting
  ## Extract distances and group info
  distances <- bd$distances
  group <- factor(bd$group)
  
  groups <- c("Control_0", "Control_7", "Control_14", "Control_21", 
              "Control_28", "Control_35", "Control_37", "Control_39",
              "Exposed_0", "Exposed_7", "Exposed_14", "Exposed_21", 
              "Exposed_28", "Exposed_35", "Exposed_37", "Exposed_39")
  
  ## P-values matrix
  ### Create empty matrices 
  mymat <- matrix(NA, ncol=16 , nrow=16)
  colnames(mymat) <- rownames(mymat) <- groups
  
  ### Fill in matrix
  mymat[lower.tri(mymat, diag = FALSE)] <- as.matrix(as.numeric(unlist(lapply(Mod_post[2:121], function(x) x['Pr(>F)'][[1]][1]))))
  mymat <- t(mymat)
  mymat[lower.tri(mymat, diag = FALSE)] <- as.matrix(as.numeric(tukey$group[,4]))
  
  # Set breaks
  mybreaks <- c(10^(seq(log10(min(mymat, na.rm=T)), log10(0.05), length.out=3)), 1)
  mybreaks <- trunc(mybreaks*10^3)/10^3 # Cut off decimals to 3 spaces
  
  ### Plot matrix as heatmap
  plot <- plotBDperm(pvalue.matrix = mymat, mybreaks = mybreaks, groups = groups, font.size=5) +
    xlab(expression(paste(beta,"-dispersion test"))) + ylab("Permanova test") + ggtitle(paste0(tissue, ": ", distance_method))
  
  # Results
  results <- list(
    Tissue = tissue,
    Method = distance_method,
    BD = bd,
    BD_Tukey = tukey,
    BD_ANOVA = ano,
    BD_Permu = permu,
    Plot = plot
  )
  return(results)
}



####---Plot distance ordination---####
#------------------------------------#

plotOrd <- function(data, dist_method, dist_gr, 
                    color_var, shape_var, ellipse_var, plot_title, ...) {
  
  # Calculate ordination
  pcoa <- ordinate(data, "PCoA", distance = dist_method)
  
  # Extract variable data
  color_var_data <- sample_data(data)[[color_var]]
  shape_var_data <- sample_data(data)[[shape_var]]
  ellipse_var_data <- sample_data(data)[[ellipse_var]]
  
  # Plotting
  p1 <- plot_ordination(dist_gr, pcoa, color = color_var_data, shape = shape_var_data, axes = c(1,2))
  p2 <- plot_ordination(dist_gr, pcoa, color = color_var_data, shape = shape_var_data, axes = c(1,3))
  p3 <- plot_ordination(dist_gr, pcoa, color = color_var_data, shape = shape_var_data, axes = c(3,2))
  
  ## Convert to ggplot for customization
  g1 <- ggplot(p1$data, p1$mapping) + 
    geom_point(aes(color = color_var_data, shape = shape_var_data), size = 3, alpha = 0.5) +
    theme(legend.key.width = unit(2, 'cm'), axis.title = element_text(face = "bold")) +
    ggtitle("") +
    xlab(p1$labels$x) +
    ylab(p1$labels$y) +
    labs(color = paste0(color_var), shape = paste0(shape_var)) +
    stat_ellipse(aes(group = ellipse_var_data, color = ellipse_var_data), level = 0.95)
  g2 <- ggplot(p2$data, p2$mapping) + 
    geom_point(aes(color = color_var_data, shape = shape_var_data), size = 3, alpha = 0.5) +
    theme(legend.key.width = unit(2, 'cm'), plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(face = "bold")) +
    ggtitle(paste(plot_title, sep = "")) +
    xlab(p2$labels$x) +
    ylab(p2$labels$y) +
    stat_ellipse(aes(group = ellipse_var_data, color = ellipse_var_data), level = 0.95)
  g3 <- ggplot(p3$data, p3$mapping) + 
    geom_point(aes(color = color_var_data, shape = shape_var_data), size = 3, alpha = 0.5) +
    theme(legend.key.width = unit(2, 'cm'), axis.title = element_text(face = "bold")) +
    ggtitle("") +
    xlab(p3$labels$x) +
    ylab(p3$labels$y) +
    stat_ellipse(aes(group = ellipse_var_data, color = ellipse_var_data), level = 0.95)
  
  ## arrange plots
  plots <- ggarrange(g1, g2, g3, ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")
  
  #Print for pdf file
  print(plots) 
}

