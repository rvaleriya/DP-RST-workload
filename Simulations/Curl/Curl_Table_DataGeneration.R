library(gt)

setwd("~/Desktop/DP-RST-workload")

#-------------------------------------------------------------------------------

mu_labels <- c("μ₁", "μ₂", "μ₃", "μ₄", "μ₅", "μ₆", "μ₇", "μ₈", "μ₉", "μ₁₀")

# Hardcoded values
table_data <- data.frame(
  Feature = paste0("PC", 1:10),
  C1 = c(0.430, -0.375, 0.803, -0.351, -1.489, -1.295, 0.123, 1.159, 0.349, 0.149),
  C2 = c(0.335, 0.187, -1.277, 0.596, 1.240, 0.768, -1.050, -1.211, -1.253, -0.256),
  C3 = c(-0.561, -0.729, -0.312, -1.459, -1.422, -1.241, 0.077, 0.555, 1.926, 0.849),
  C4 = c(0.897, 0.380, 0.935, -0.160, -1.355, -0.652, -1.188, 0.761, 1.270, -1.743),
  C5 = c(-0.749, -1.062, 0.060, 1.535, -0.203, -0.705, -0.606, -0.775, 1.620, 0.913),
  C6 = c(-0.770, 0.734, 1.417, 1.063, -1.373, -0.970, 1.616, 1.969, -1.779, 1.587),
  C7 = c(-0.430, -1.479, 0.631, 1.097, -0.918, -0.230, -1.795, 1.555, 1.631, -0.804),
  C8 = c(-1.029, 1.397, -0.845, -0.961, -1.711, 0.713, 0.935, -1.869, 0.461, -1.540),
  C9 = c(1.540, -1.572, -1.959, -1.951, 1.779, 1.854, 1.320, -1.995, -0.610, 1.877),
  C10 = c(-1.703, 0.915, -0.465, 1.461, 0.073, 1.613, 0.978, 0.283, -0.492, -1.984),
  Scale = c(0.891, 0.581, 0.835, 2.189, 1.905, 1.581, 2.481, 1.159, 1.088, 2.592),
  Slant = c(3.345, -8.479, 3.540, -4.110, -8.315, -3.754, -8.654, 7.299, 7.435, 9.930)
)

#-------------------------------------------------------------------------------

gt_table <- table_data %>%
  gt() %>%
  fmt_number(everything(), decimals = 3) %>%
  tab_spanner(
    label = md("**Cluster Means**"),
    columns = starts_with("C")
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = "Feature")
  ) %>%
  cols_label(
    Feature = md("**Feature**"),
    C1 = html("<b>μ₁</b>"),
    C2 = html("<b>μ₂</b>"),
    C3 = html("<b>μ₃</b>"),
    C4 = html("<b>μ₄</b>"),
    C5 = html("<b>μ₅</b>"),
    C6 = html("<b>μ₆</b>"),
    C7 = html("<b>μ₇</b>"),
    C8 = html("<b>μ₈</b>"),
    C9 = html("<b>μ₉</b>"),
    C10 = html("<b>μ₁₀</b>"),
    Scale = md("**Scale**"),
    Slant = md("**Slant**")
  ) %>%
  cols_width(
    everything() ~ px(80)
  ) %>%
  cols_align(
    align = "center",
    columns = c(starts_with("C"), Scale, Slant)
  ) %>%
  tab_style(
    style = cell_borders(sides = "left", color = "#C0C0C0", weight = px(2)),
    locations = list(
      cells_body(columns = "C1"),
      cells_body(columns = "Scale"),
      cells_body(columns = "Slant")
    )
  )

gt_table

#-------------------------------------------------------------------------------

gtsave(gt_table, "New_Simulations/Curl/cluster_table.png",  
       vwidth = 1600, vheight = 900, dpi = 600)

#-------------------------------------------------------------------------------

library(ggplot2)
library(cowplot)
library(magick)

img <- magick::image_read("New_Simulations/Curl/cluster_table.png")

labeled_table <- ggdraw() +
  draw_image(img, scale = 1) +
  draw_label("(B)", x = 0.01, y = 0.99, hjust = 0, vjust = 1, size = 20)

labeled_table

ggsave("New_Simulations/Curl/cluster_table_labeled_C.png", plot = labeled_table,
       width = 10, height = 5, dpi = 300)



