# Load necessary libraries
library(ggplot2)
library(sp)
library(plotly)
library(rprojroot)

# Set working directory to the root of the R project
project_root <- find_root(is_rstudio_project)
setwd(project_root)
getwd()

# Load datasets
load("./Real Data/swiss_roll_wt_muscle_finaltouches1.RData")
load("./Real Data/swiss_roll_wt_muscle_boundary.RData")
load("./Real Data/Swiss_Roll_DPM_OutputOnly.RData")

# Extract first three principal components and spatial coordinates
Y_sample <- swiss_roll_wt_muscle_finaltouches1[, 1:3]
loc <- swiss_roll_wt_muscle_finaltouches1[, 4:5]
bnd <- boundary

# Standardize coordinates and Y values
coords <- as.data.frame(apply(loc, 2, scale))
Y_std <- as.data.frame(apply(Y_sample, 2, scale))

# Standardize the boundary given the params of coordinates
bnd_scaled <- data.frame(
  x = (bnd$x - mean(loc[, 1])) / sd(loc[, 1]),
  y = (bnd$y - mean(loc[, 2])) / sd(loc[, 2])
)

# Add small random variable to the scaled boundary to avoid duplicates
set.seed(123)  # For reproducibility
bnd_scaled_r <- data.frame(
  x = bnd_scaled$x + c(0, runif(length(bnd_scaled$x) - 2, min = 1e-04, max = 5e-04), 0),
  y = bnd_scaled$y + c(0, runif(length(bnd_scaled$y) - 2, min = 1e-04, max = 5e-04), 0)
)

# Create dataframe for plotting
df_subset <- data.frame(coords, Y_std)

# Get the palette from paletteer
palette_colors <- paletteer::paletteer_c("grDevices::Blue-Red 2", n = 100)

# Convert the palette to a format usable by Plotly
palette_colors <- as.character(palette_colors)

# Create a colorscale that Plotly can use
colorscale <- list()
for (i in seq_along(palette_colors)) {
  colorscale[[i]] <- list((i - 1) / (length(palette_colors) - 1), palette_colors[i])
}

# Create the plot using plotly
p <- plot_ly() %>%
  add_trace(
    data = df_subset,
    x = ~(-x), z = ~y, y = 0,
    type = 'scatter3d', mode = 'markers',
    marker = list(size = 3, color = ~PC3, colorscale = colorscale, opacity = 1, 
                  colorbar = list(len = 0.8, title = "PCA values")),
    showlegend = FALSE  # Hide legend
  ) %>%
  add_trace(
    data = df_subset,
    x = ~(-x), z = ~y, y = 5,
    type = 'scatter3d', mode = 'markers',
    marker = list(size = 3, color = ~PC2, colorscale = colorscale, opacity = 1),
    showlegend = FALSE  # Hide legend
  ) %>%
  add_trace(
    data = df_subset,
    x = ~(-x), z = ~y, y = 8,
    type = 'scatter3d', mode = 'markers',
    marker = list(size = 3, color = ~PC1, colorscale = colorscale, opacity = 1),
    showlegend = FALSE  # Hide legend
  ) %>%
  add_trace(
    data = bnd_scaled_r,
    x = ~(-x), z = ~y, y = rep(0, nrow(bnd_scaled_r)),
    type = 'scatter3d', mode = 'lines',
    line = list(color = 'black', width = 2),
    showlegend = FALSE  # Hide legend
  ) %>%
  add_trace(
    data = bnd_scaled_r,
    x = ~(-x), z = ~y, y = rep(5, nrow(bnd_scaled_r)),
    type = 'scatter3d', mode = 'lines',
    line = list(color = 'black', width = 2),
    showlegend = FALSE  # Hide legend
  ) %>%
  add_trace(
    data = bnd_scaled_r,
    x = ~(-x), z = ~y, y = rep(8, nrow(bnd_scaled_r)),
    type = 'scatter3d', mode = 'lines',
    line = list(color = 'black', width = 2),
    showlegend = FALSE  # Hide legend
  ) %>%
  layout(scene = list(
    xaxis = list(title = '', showticklabels = FALSE, showbackground = FALSE, zeroline = FALSE, showgrid = FALSE),
    yaxis = list(title = '', showticklabels = FALSE, showbackground = FALSE, zeroline = FALSE, showgrid = FALSE),
    zaxis = list(title = '', showticklabels = FALSE, showbackground = FALSE, zeroline = FALSE, showgrid = FALSE),
    bgcolor = 'rgba(0,0,0,0)'
  ),
  paper_bgcolor = 'rgba(0,0,0,0)',
  plot_bgcolor = 'rgba(0,0,0,0)'
  )

p
