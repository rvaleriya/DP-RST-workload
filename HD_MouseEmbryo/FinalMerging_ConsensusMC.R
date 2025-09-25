
set.seed(42)

load("~/Desktop/DP-RST-workload/HD_Mouse_Embryo/Processed_Data/MouseEmbryo_HD_merged_data.RData")

local_clusters_and_anchor_idx <- readRDS("~/Desktop/DP-RST-workload/HD_Mouse_Embryo/local_clusters_and_anchor_idx.rds")

local_clusters <- local_clusters_and_anchor_idx[["local_clusters"]]
anchor_idx <- local_clusters_and_anchor_idx[["anchor_idx"]]

calculate_anchor_distance <- function(cluster1_indices, cluster2_indices, anchor_idx) {
  anchors_in_cluster1 <- intersect(cluster1_indices, anchor_idx)
  anchors_in_cluster2 <- intersect(cluster2_indices, anchor_idx)
  
  common_anchors_count <- length(intersect(anchors_in_cluster1, anchors_in_cluster2))
  
  # D: Number of anchor points present in one cluster but not the other (symmetric difference)
  # diff1 = anchors_in_cluster1 that are NOT in anchors_in_cluster2
  # diff2 = anchors_in_cluster2 that are NOT in anchors_in_cluster1
  different_anchors_count <- length(setdiff(anchors_in_cluster1, anchors_in_cluster2)) +
    length(setdiff(anchors_in_cluster2, anchors_in_cluster1))
  
  # Handle edge cases as per paper: if C=0 and D=0 (e.g. both clusters have no anchors) distance is 1.
  if (common_anchors_count == 0 && different_anchors_count == 0) {
    return(1.0) 
  }
  # If C > 0 and D = 0 (identical non-empty anchor sets), distance is 0.
  if (different_anchors_count == 0 && common_anchors_count > 0) {
    return(0.0)
  }
  # If C = 0 and D > 0 (disjoint non-empty anchor sets, or one empty other non-empty) distance is 1
  if (common_anchors_count == 0 && different_anchors_count > 0) {
    return(1.0)
  }
  
  distance <- different_anchors_count / (common_anchors_count + different_anchors_count)
  return(distance)
}

merge_clusters_by_anchors <- function(local_clusters_list, anchor_indices, epsilon_threshold) {
  if (!is.list(local_clusters_list) || length(local_clusters_list) == 0) {
    return(local_clusters_list)
  }
  
  # Diagnostic: Check clusters with no anchors
  num_no_anchors <- 0
  num_total_anchors_in_local_clusters <- 0
  min_anchors_in_cluster <- Inf
  max_anchors_in_cluster <- 0
  
  for (cl in local_clusters_list) {
    anchors_present <- intersect(cl, anchor_indices)
    if (length(anchors_present) == 0) {
      num_no_anchors <- num_no_anchors + 1
    }
    num_total_anchors_in_local_clusters <- num_total_anchors_in_local_clusters + length(anchors_present)
    if (length(local_clusters_list) > 0) { # only update if list is not empty
      min_anchors_in_cluster <- min(min_anchors_in_cluster, length(anchors_present))
      max_anchors_in_cluster <- max(max_anchors_in_cluster, length(anchors_present))
    } else {
      min_anchors_in_cluster <- 0 # if list is empty, no anchors
    }
    
  }
  cat("Initial number of local clusters:", length(local_clusters_list), "\n")
  cat("Number of local clusters with NO anchor points:", num_no_anchors, "\n")
  cat("Average anchors per local cluster (if any):", 
      ifelse(length(local_clusters_list) - num_no_anchors > 0, 
             num_total_anchors_in_local_clusters / (length(local_clusters_list) - num_no_anchors), 
             0), "\n")
  cat("Min anchors in a local cluster:", ifelse(is.infinite(min_anchors_in_cluster), 0, min_anchors_in_cluster), "\n")
  cat("Max anchors in a local cluster:", max_anchors_in_cluster, "\n")
  cat("Total unique anchor points provided:", length(unique(anchor_indices)), "\n")
  
  if (num_no_anchors == length(local_clusters_list) && length(local_clusters_list) > 0) {
    cat("Warning: All local clusters have no anchor points. Cannot merge based on anchors.\n")
    return(local_clusters_list)
  }
  
  global_clusters <- local_clusters_list
  iteration_count <- 0
  cat("Starting merge process with epsilon =", epsilon_threshold, "...\n")
  
  while (TRUE) {
    iteration_count <- iteration_count + 1
    cat("Iteration:", iteration_count, ", Current number of clusters:", length(global_clusters), "\n")
    
    if (length(global_clusters) <= 1) {
      cat("Only one or zero clusters left. Stopping.\n")
      break
    }
    
    # Randomly permute the order of clusters (as per paper's Algorithm 1)
    if (length(global_clusters) > 1) {
      global_clusters <- global_clusters[sample(length(global_clusters))]
    }
    
    merged_in_this_pass <- FALSE
    temp_global_clusters <- list() 
    
    # Keep track of which clusters (by index in current 'global_clusters' list) have been processed in this pass
    processed_indices_this_pass <- logical(length(global_clusters)) 
    
    i <- 1
    while(i <= length(global_clusters)){
      if(processed_indices_this_pass[i]){ # Skip if already merged/processed in this pass
        i <- i + 1
        next
      }
      
      clusterA_content <- global_clusters[[i]]
      best_merge_candidate_idx_in_current_list <- -1
      min_distance_for_A <- Inf
      
      # Search for a merge partner for clusterA
      if (i < length(global_clusters)) {
        for (j in (i + 1):length(global_clusters)) {
          if(processed_indices_this_pass[j]){ # Skip if partner already merged/processed
            next
          }
          
          clusterB_content <- global_clusters[[j]]
          distance <- calculate_anchor_distance(clusterA_content, clusterB_content, anchor_indices)
          
          if (distance < epsilon_threshold && distance < min_distance_for_A) {
            min_distance_for_A <- distance
            best_merge_candidate_idx_in_current_list <- j
          }
        }
      }
      
      if (best_merge_candidate_idx_in_current_list != -1) {
        cluster_to_merge_with_content <- global_clusters[[best_merge_candidate_idx_in_current_list]]
        # Perform the merge: unique points from both clusters
        merged_cluster <- unique(c(clusterA_content, cluster_to_merge_with_content))
        temp_global_clusters <- c(temp_global_clusters, list(merged_cluster))
        
        # cat("  Merged cluster (orig_idx unknown after shuffle, current list pos", i, 
        #     ") with (current list pos", best_merge_candidate_idx_in_current_list, 
        #     ") (distance:", sprintf("%.3f", min_distance_for_A), 
        #     ") New size:", length(merged_cluster), "\n")
        
        processed_indices_this_pass[i] <- TRUE
        processed_indices_this_pass[best_merge_candidate_idx_in_current_list] <- TRUE
        merged_in_this_pass <- TRUE
      } else {
        # If clusterA was not merged, add it to the list for the next pass
        temp_global_clusters <- c(temp_global_clusters, list(clusterA_content))
        processed_indices_this_pass[i] <- TRUE # Mark as processed (i.e., kept as is for this pass)
      }
      i <- i + 1 
    } # End of while(i <= length(global_clusters))
    
    global_clusters <- temp_global_clusters 
    
    if (!merged_in_this_pass) {
      cat("No merges in this pass. Convergence reached.\n")
      break
    }
    if (iteration_count > 200) { # Increased safety break
      cat("Warning: Exceeded 200 iterations. Breaking loop.\n")
      break
    }
  } # End of while(TRUE)
  
  cat("Merge process finished. Final number of clusters:", length(global_clusters), "\n")
  return(global_clusters)
}

epsilon <- 0.2 # Example, adjust based on results
final_clusters_high_eps <- merge_clusters_by_anchors(local_clusters, 
                                                     anchor_idx, 
                                                     epsilon)
print(paste("Number of final clusters:", length(final_clusters_high_eps)))

#-------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# Optimized Cluster Visualization Code
# --------------------------------------------------------------------------------

# Load required packages upfront
required_packages <- c("ggplot2", "dplyr", "RColorBrewer", "viridisLite", "scales")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# --------------------------------------------------------------------------------
# Helper Functions
# --------------------------------------------------------------------------------

#' Get a color palette based on the number of clusters
#' 
#' @param num_clusters Number of clusters to create colors for
#' @return A vector of colors
get_cluster_palette <- function(num_clusters) {
  if (num_clusters <= 0) return(character(0))
  
  if (num_clusters <= 12) {
    palette_name <- if (num_clusters < 3) "Set1" else "Paired"
    palette_values <- brewer.pal(max(3, num_clusters), palette_name)
    if (num_clusters < length(palette_values)) {
      palette_values <- palette_values[1:num_clusters]
    }
  } else {
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
    
    if (num_clusters <= length(col_vector)) {
      set.seed(123) # for reproducibility
      palette_values <- sample(col_vector, num_clusters)
    } else {
      palette_values <- viridis(num_clusters, option = "D")
    }
  }
  
  return(palette_values)
}

#' Create a data frame from cluster indices
#' 
#' @param cluster_list List of clusters with point indices
#' @param data_df Original data frame with coordinates
#' @param x_col Name of x coordinate column
#' @param y_col Name of y coordinate column
#' @param selected_clusters Optional vector of cluster indices to include
#' @return Data frame with points and their cluster assignments
prepare_cluster_data <- function(cluster_list, data_df, x_col = "x", y_col = "y", 
                                 selected_clusters = NULL) {
  
  # Validate coordinate columns
  if (!all(c(x_col, y_col) %in% names(data_df))) {
    stop(paste("Coordinate columns '", x_col, "' or '", y_col, "' not found in data."))
  }
  
  # If no specific clusters are selected, use all
  if (is.null(selected_clusters)) {
    selected_clusters <- 1:length(cluster_list)
  }
  
  # Create a list to store data frames for each cluster
  cluster_dfs <- list()
  
  for (i in seq_along(selected_clusters)) {
    cluster_id <- selected_clusters[i]
    point_indices <- cluster_list[[cluster_id]]
    
    # Ensure indices are valid
    valid_indices <- point_indices[point_indices > 0 & point_indices <= nrow(data_df)]
    
    if (length(valid_indices) > 0) {
      # Get relevant data for these points
      cluster_points <- data_df[valid_indices, c(x_col, y_col), drop = FALSE]
      cluster_points$cluster_id <- as.factor(i)  # Use sequential IDs for plotting
      cluster_points$original_cluster_id <- as.factor(cluster_id)  # Keep original ID for reference
      
      cluster_dfs[[i]] <- cluster_points
    }
  }
  
  # Combine all clusters into one data frame
  if (length(cluster_dfs) > 0) {
    return(do.call(rbind, cluster_dfs))
  } else {
    return(NULL)
  }
}

#' Create a cluster visualization plot
#' 
#' @param data_df Data frame with points and cluster assignments
#' @param x_col Name of x coordinate column
#' @param y_col Name of y coordinate column
#' @param cluster_col Name of cluster ID column
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param colors Optional color palette
#' @param point_size Size of points
#' @param point_alpha Alpha transparency of points
#' @param facet Whether to create a faceted plot
#' @return ggplot object
plot_cluster_visualization <- function(data_df, x_col = "x", y_col = "y", 
                                       cluster_col = "cluster_id", 
                                       title = "Cluster Visualization",
                                       subtitle = NULL,
                                       colors = NULL,
                                       point_size = 1.5,
                                       point_alpha = 0.7,
                                       facet = FALSE) {
  
  if (is.null(data_df) || nrow(data_df) == 0) {
    message("No data to plot.")
    return(NULL)
  }
  
  num_clusters <- length(unique(data_df[[cluster_col]]))
  
  # if (is.null(subtitle)) {
  #   subtitle <- paste("Total points:", nrow(data_df), "in", num_clusters, "clusters")
  # }
  
  # Base plot
  p <- ggplot(data_df, aes_string(x = x_col, y = y_col, color = cluster_col)) +
    geom_point(size = point_size, alpha = point_alpha) +
    theme_minimal(base_size = 12) +
    labs(title = title,
         subtitle = subtitle,
         x = gsub("_", " ", tools::toTitleCase(x_col)),
         y = gsub("_", " ", tools::toTitleCase(y_col)),
         color = "Cluster ID") +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 16))
  
  # Add colors if provided, otherwise generate them
  if (is.null(colors) && num_clusters > 0) {
    colors <- get_cluster_palette(num_clusters)
  }
  
  if (length(colors) > 0) {
    p <- p + scale_color_manual(values = colors, na.value = "grey80")
  } else if (num_clusters > 0) {
    p <- p + scale_color_viridis_d(option = "turbo", na.value = "grey80")
  }
  
  # Create faceted plot if requested
  if (facet && num_clusters > 1 && num_clusters <= 25) {
    p <- p + 
      facet_wrap(as.formula(paste("~", cluster_col)), 
                 ncol = floor(sqrt(num_clusters))) +
      theme(legend.position = "none")
  }
  
  return(p)
}

# --------------------------------------------------------------------------------
# Main Analysis Functions
# --------------------------------------------------------------------------------

#' Analyze and plot cluster sizes
#' 
#' @param clusters List of clusters with point indices
#' @return Data frame with cluster size information
analyze_cluster_sizes <- function(clusters) {
  # Calculate cluster sizes
  cluster_sizes <- sapply(clusters, length)
  
  # Create summary data frame
  cluster_summary <- data.frame(
    cluster_id = 1:length(cluster_sizes),
    size = cluster_sizes
  )
  
  # Print summary statistics
  cat("\nSizes of each cluster:\n")
  for (i in 1:length(cluster_sizes)) {
    cat("Cluster", i, ":", cluster_sizes[i], "observations\n")
  }
  
  cat("\nSummary of cluster sizes:\n")
  print(summary(cluster_sizes))
  
  # Create boxplot of cluster sizes
  p <- ggplot(cluster_summary, aes(y = size)) +
    geom_boxplot(fill = "skyblue") +
    theme_minimal() +
    labs(title = "Boxplot of Cluster Sizes",
         y = "Number of Observations per Cluster") +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )
  
  print(p)
  
  return(cluster_summary)
}

#' Plot the spatial distribution of all clusters
#' 
#' @param clusters List of clusters with point indices
#' @param data_df Original data frame with coordinates
#' @param x_col Name of x coordinate column
#' @param y_col Name of y coordinate column
#' @param epsilon Value of epsilon parameter (optional)
#' @return ggplot object
plot_all_clusters <- function(clusters, data_df, x_col = "x", y_col = "y", epsilon = NULL) {
  # Prepare data
  plot_data <- prepare_cluster_data(clusters, data_df, x_col, y_col)
  
  # Create title
  title <- "Spatial Distribution of All Clusters"
  if (!is.null(epsilon)) {
    title <- paste(title, "(Epsilon =", epsilon, ")")
  }
  
  # Create and return plot
  plot_cluster_visualization(
    plot_data, 
    x_col = x_col, 
    y_col = y_col,
    title = title
  )
}

#' Plot only the top N largest clusters
#' 
#' @param clusters List of clusters with point indices
#' @param data_df Original data frame with coordinates
#' @param x_col Name of x coordinate column
#' @param y_col Name of y coordinate column
#' @param n Number of top clusters to include
#' @param epsilon Value of epsilon parameter (optional)
#' @param create_facet_plot Whether to create a faceted plot in addition to main plot
#' @return List with plot objects
plot_top_clusters <- function(clusters, data_df, x_col = "x", y_col = "y", 
                              n = 10, epsilon = NULL, create_facet_plot = TRUE) {
  
  # Calculate cluster sizes
  cluster_sizes <- sapply(clusters, length)
  cluster_summary <- data.frame(
    cluster_id = 1:length(cluster_sizes),
    size = cluster_sizes
  )
  
  # Select top N clusters
  n <- min(n, nrow(cluster_summary))
  top_clusters <- cluster_summary %>%
    arrange(desc(size)) %>%
    slice_head(n = n) %>%
    pull(cluster_id)
  
  # Report selected clusters
  cat("Selected Top", length(top_clusters), "largest clusters. IDs:", 
      paste(top_clusters, collapse = ", "), "\n")
  cat("Their sizes are:", 
      paste(cluster_summary$size[top_clusters], collapse = ", "), "\n")
  
  # Prepare data for selected clusters
  plot_data <- prepare_cluster_data(clusters, data_df, x_col, y_col, top_clusters)
  
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    message("No valid data to plot after selecting top clusters.")
    return(NULL)
  }
  
  # Get color palette for all plots
  num_clusters <- length(top_clusters)
  colors <- get_cluster_palette(num_clusters)
  
  # Create title
  title <- paste("Top", num_clusters, "Largest Clusters")
  # if (!is.null(epsilon)) {
  #   title <- paste(title, "(Epsilon =", epsilon, ")")
  # }
  
  # Create main plot
  main_plot <- plot_cluster_visualization(
    plot_data, 
    x_col = x_col, 
    y_col = y_col,
    title = title,
    colors = colors
  )
  
  result <- list(main_plot = main_plot)
  
  # Create faceted plot if requested
  if (create_facet_plot && num_clusters > 1 && num_clusters <= 25) {
    facet_plot <- plot_cluster_visualization(
      plot_data, 
      x_col = x_col, 
      y_col = y_col,
      title = paste("Faceted View:", title),
      colors = colors,
      point_size = 0.8,
      point_alpha = 0.6,
      facet = TRUE
    )
    
    print(facet_plot)
    result$facet_plot <- facet_plot
  }
  
  print(main_plot)
  return(result)
}

# --------------------------------------------------------------------------------
# Example Usage
# --------------------------------------------------------------------------------

# Analyze cluster sizes
cluster_summary <- analyze_cluster_sizes(final_clusters_high_eps)

# Plot all clusters
all_clusters_plot <- plot_all_clusters(
  final_clusters_high_eps,
  merged_data,
  x_col = "x",
  y_col = "y",
  epsilon = epsilon
)
all_clusters_plot

# Plot top 10 largest clusters
top_clusters_plots <- plot_top_clusters(
  final_clusters_high_eps,
  merged_data,
  x_col = "x",
  y_col = "y",
  n = 3,
  epsilon = epsilon,
  create_facet_plot = TRUE
)
top_clusters_plots





