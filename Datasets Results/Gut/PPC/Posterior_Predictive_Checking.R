# ===================================================================
# Load Libraries & Data
# ===================================================================
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(cluster) 
library(grid)

setwd("~/Desktop/DP-RST-workload/Datasets Results/Gut")

set.seed(123)

# Load data and results
gut_df_wt_muscle <- readRDS("gut_df_wt_muscle.rds")

# ===================================================================
# MAIN PPC ANALYSIS FUNCTION 
# ===================================================================
run_simplified_ppc <- function(gut_df_wt_muscle, output, n_pcs, output_suffix) {
  
  cat("\n===== Running PPC Analysis for PC Means -", n_pcs, "PCs =====\n\n")
  
  # Extract the observed data (Principal Components)
  Y_obs <- as.matrix(gut_df_wt_muscle[, 1:n_pcs])
  Y_obs <- apply(Y_obs, 2, scale)
  n <- nrow(Y_obs)
  p <- ncol(Y_obs)
  S <- length(output$mu_out) # Number of MCMC samples
  
  # ===================================================================
  # Generate All Replicated Datasets (FROM TEAMS/REFINED PARTITION)
  # ===================================================================
  cat("Generating", S, "replicated datasets from TEAMS partition...\n")
  
  replicated_datasets_list <- lapply(1:S, function(s) {
    if (s %% 500 == 0) cat("  Generating dataset", s, "of", S, "...\n")
    
    # Get the TEAMS parameters
    mu_teams_s <- output$mu_teams_out[[s]]  # Teams means
    Sigma_s <- output$sigmasq_y_out[[s]]    # Covariance (same for all teams)
    
    # Get the REFINED partition (teams assignment)
    spatial_partition_s <- output$cluster_out[s, ]
    refined_mapping_s <- output$teams_out[[s]]
    refined_partition_s <- refined_mapping_s[spatial_partition_s]
    
    Y_rep_s <- matrix(NA, nrow = n, ncol = p)
    
    # Generate data for each TEAM (not spatial cluster)
    for (team in unique(refined_partition_s)) {
      indices <- which(refined_partition_s == team)
      Y_rep_s[indices, ] <- rmvnorm(length(indices), 
                                    mean = mu_teams_s[team, ], 
                                    sigma = Sigma_s)
    }
    return(Y_rep_s)
  })
  
  cat("Done generating all replicated datasets.\n\n")
  
  # ===================================================================
  # Calculate PC Means 
  # ===================================================================
  cat("Calculating PC means...\n")
  
  # Function for PC means
  compute_means <- function(Y) {
    p <- ncol(Y)
    means <- sapply(1:p, function(j) mean(Y[, j]))
    names(means) <- paste0("mean_PC", 1:p)
    return(means)
  }
  
  # Calculate means for observed data
  obs_means <- compute_means(Y_obs)
  
  # Calculate means for replicated data
  rep_means_df <- as.data.frame(t(sapply(replicated_datasets_list, compute_means)))
  
  cat("Done calculating statistics.\n\n")
  
  # ===================================================================
  # Calculate PPP-values for Means
  # ===================================================================
  
  ppp_means <- sapply(names(obs_means), function(stat) {
    ppp <- mean(rep_means_df[[stat]] >= obs_means[stat], na.rm = TRUE)
    return(min(ppp, 1-ppp))  # Two-tailed extreme
  })
  
  # Print summary
  cat("PPP VALUES FOR PC MEANS (values < 0.025 indicate potential misfit):\n")
  for(i in 1:length(ppp_means)) {
    flag <- if(ppp_means[i] < 0.025) " ⚠️" else " ✓"
    cat(sprintf("  PC%d mean: %.3f%s\n", i, ppp_means[i], flag))
  }
  
  # ===================================================================
  # Create Plots for Means
  # ===================================================================
  cat("\nGenerating means plots...\n")
  
  means_plots <- lapply(1:length(obs_means), function(i) {
    stat <- names(obs_means)[i]
    ppp_val <- ppp_means[stat]
    
    pc_num <- i
    
    p <- ggplot(rep_means_df, aes_string(x = stat)) +
      geom_histogram(bins = 40, fill = "skyblue", alpha = 0.8, color = "white") +
      geom_vline(xintercept = obs_means[stat], 
                 color = "red", linetype = "dashed", linewidth = 1) +
      labs(title = paste0("PC", pc_num, " Mean"),
           subtitle = sprintf("PPP = %.3f", ppp_val),
           x = "Mean Value", 
           y = "Frequency") +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 11, color = "gray40"))
    
    # Add red background if misfit
    if (ppp_val < 0.025) {
      p <- p + theme(panel.background = element_rect(fill = "#ffe6e6", color = NA))
    }
    p
  })
  
  # ===================================================================
  # Save Plots Based on Number of PCs
  # ===================================================================
  
  if (n_pcs == 3) {
    # For 3 PCs: 1 row, 3 columns
    pdf(paste0("PPC/PPC_Means_", output_suffix, ".pdf"), width = 12, height = 4)
    grid.arrange(grobs = means_plots, ncol = 3, nrow = 1)
    dev.off()
    
  } else if (n_pcs == 10) {
    # For 10 PCs: 2 rows, 5 columns
    pdf(paste0("PPC/PPC_Means_", output_suffix, ".pdf"), width = 15, height = 6)
    grid.arrange(grobs = means_plots, ncol = 5, nrow = 2)
    dev.off()
  }
  
  cat(sprintf("Plot saved: PPC/PPC_Means_%s.pdf\n", output_suffix))
  
  # Return results
  return(list(ppp_means = ppp_means))
}

# ===================================================================
# RUN ANALYSIS FOR 3 PCs
# ===================================================================
cat("\n################################################\n")
cat("ANALYSIS 1: 3 Principal Components\n")
cat("################################################\n")

# Create PPC directory if it doesn't exist
if (!dir.exists("PPC")) dir.create("PPC")

# Load 3 PC results
load("Results/Gut_DP.RST_FromNewBastPT_p3_Version2_OutputOnly.RData")
output_3pc <- Gut_DP.RST_FromNewBastPT_Version2

# Run analysis
results_3pc <- run_simplified_ppc(
  gut_df_wt_muscle = gut_df_wt_muscle,
  output = output_3pc,
  n_pcs = 3,
  output_suffix = "3PC"
)

# ===================================================================
# RUN ANALYSIS FOR 10 PCs
# ===================================================================
cat("\n################################################\n")
cat("ANALYSIS 2: 10 Principal Components\n")
cat("################################################\n")

# Load 10 PC results
load("Results/Gut_DP.RST_FromNewBastPT_p10_Version2_OutputOnly.RData")
output_10pc <- Gut_DP.RST_FromNewBastPT_p10_Version2

# Run analysis
results_10pc <- run_simplified_ppc(
  gut_df_wt_muscle = gut_df_wt_muscle,
  output = output_10pc,
  n_pcs = 10,
  output_suffix = "10PC"
)

# ===================================================================
# FINAL SUMMARY
# ===================================================================
cat("\n################################################\n")
cat("PPC ANALYSIS COMPLETE\n")
cat("################################################\n\n")

# Count issues
issues_3pc <- sum(results_3pc$ppp_means < 0.025)
issues_10pc <- sum(results_10pc$ppp_means < 0.025)

cat("Summary of PC Means:\n")
cat(sprintf("  3 PC Model:  %d PC means with potential misfit\n", issues_3pc))
cat(sprintf("  10 PC Model: %d PC means with potential misfit\n", issues_10pc))

if (issues_3pc < issues_10pc) {
  cat("\n→ The 3 PC model shows better fit for PC means.\n")
} else if (issues_10pc < issues_3pc) {
  cat("\n→ The 10 PC model shows better fit for PC means.\n")
} else {
  cat("\n→ Both models show similar fit for PC means.\n")
}

cat("\nFiles generated:\n")
cat("  • PPC/PPC_Means_3PC.pdf (3 plots in 1 row)\n")
cat("  • PPC/PPC_Means_10PC.pdf (10 plots in 2×5 grid)\n")