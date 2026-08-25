# ==============================================================================
# R/utils.R: Helper Functions for Posterior Summaries and Visualization
# ==============================================================================

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggdist))

# 1. Summary Statistics Helpers ------------------------------------------------

summarise_or <- function(df, value_col, group_cols) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      mean = mean(.data[[value_col]]),
      median = median(.data[[value_col]]),
      lower = quantile(.data[[value_col]], 0.025),
      upper = quantile(.data[[value_col]], 0.975),
      .groups = "drop"
    )
}

summarise_plot_stats <- function(draws_df) {
  draws_df %>%
    group_by(condition, eq) %>%
    summarise(
      mean = mean(or),
      median = median(or),
      sd = sd(or),
      q2_5 = quantile(or, 0.025),
      q97_5 = quantile(or, 0.975),
      .groups = "drop"
    )
}

# 2. Table Labeling Helpers ----------------------------------------------------

factor_label_for_table <- function(x) {
  label_map <- c(
    "Enjoyment" = "Enjoyment\n(vs. Null)",
    "Distraction" = "Distraction\n(vs. Null)",
    "Relaxation" = "Relaxation\n(vs. Null)",
    "Focus-Motivation" = "Focus-Motivation\n(vs. Null)",
    "DIESEL" = "Diesel\n(vs. EV)",
    "Content" = "Music Excerpt\n",
    "Individual" = "Participant\n"
  )
  recode(x, !!!label_map, .default = x) %>% gsub("\\s*\\n\\s*", " ", x = .)
}

# 3. Reusable ggplot2 Halfeye Panel Generator ----------------------------------

#' Generate a standardized halfeye coefficient distribution plot.
#'
#' @param draws_data A data frame containing posterior draws.
#' @param summary_data A data frame containing point estimates and credible intervals.
#' @param x_levels Character vector defining the factor order on the x-axis.
#' @param x_labels Named character vector for x-axis display labels.
#' @param y_limits Numeric vector of length 2 defining the y-axis Cartesian bounds.
#' @param y_breaks Numeric vector defining y-axis tick locations.
#' @param y_labels_fn Function for formatting y-axis text.
#' @param y_title String for the y-axis title.
#' @param panel_title String for the plot title.
#' @param color_map Named character vector mapping EQ profiles to hex colors.
#' @param reference_line Numeric value for the dashed null-effect line (e.g., 1 for OR, 0 for linear).
#' @param is_random_effect Boolean. If TRUE, plots the median instead of the mean.
#' @param show_y_axis_title Boolean.
#' @param show_x_axis_labels Boolean.
#' @param show_legend Boolean.
#' @return A ggplot object.

create_effect_panel <- function(
    draws_data,
    summary_data,
    x_levels,
    x_labels,
    y_limits,
    y_breaks,
    y_labels_fn,
    y_title,
    panel_title,
    color_map,
    reference_line = 1,
    is_random_effect = FALSE,
    show_y_axis_title = TRUE,
    show_x_axis_labels = TRUE,
    show_legend = FALSE
) {
  
  # Prepare summary values
  summary_data_plot <- summary_data
  if (is_random_effect) {
    summary_data_plot$point_val <- summary_data_plot$median
  } else {
    summary_data_plot$point_val <- summary_data_plot$mean
  }
  
  # Safe check for NA or NULL panel titles
  has_title <- !is.null(panel_title) && !is.na(panel_title) && panel_title != ""
  
  is_volume_model <- length(unique(draws_data$eq)) == 1 && unique(draws_data$eq)[1] %in% c("Volume", "Volume Level")
  dodge_pos <- if (is_volume_model) "identity" else position_dodge(width = 0.8)
  
  slab_w <- if (is_volume_model) 0.5 else 0.8
  
  # Force Volume plots to draw shorter densities to match the EQ visuals
  slab_scale_val <- if (is_volume_model) 0.35 else 0.9 
  
  slab_alpha_val <- if (is_random_effect) 0.8 else if (reference_line == 0) 0.5 else 0.8
  
  p <- ggplot() +
    ggdist::stat_halfeye(
      data = draws_data,
      aes(x = factor(condition, levels = x_levels), y = or, fill = eq, color = eq),
      position = dodge_pos,
      width = slab_w,
      scale = slab_scale_val,
      alpha = slab_alpha_val,
      slab_alpha = slab_alpha_val,
      interval_alpha = 0,
      point_interval = NULL,
      justification = 0
    ) +
    geom_errorbar(
      data = summary_data_plot,
      aes(x = factor(condition, levels = x_levels), ymin = lower, ymax = upper, color = eq),
      position = dodge_pos,
      width = 0.2,
      linewidth = 0.8
    ) +
    geom_point(
      data = summary_data_plot,
      aes(x = factor(condition, levels = x_levels), y = point_val, color = eq),
      position = dodge_pos,
      size = 2.5
    ) +
    geom_hline(yintercept = reference_line, linetype = "dashed", linewidth = 1.2, color = "grey40") +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(breaks = y_breaks, labels = y_labels_fn, minor_breaks = NULL) +
    coord_cartesian(ylim = y_limits) +
    scale_color_manual(values = color_map, drop = FALSE) +
    scale_fill_manual(values = color_map, drop = FALSE) +
    labs(x = NULL, y = if (show_y_axis_title) y_title else NULL)
  
  p <- p + theme_minimal(base_size = 16) +
    theme(
      axis.title.y = if (show_y_axis_title) element_text(size = 18, color = "grey40", margin = margin(r = 6, l = 4)) else element_blank(),
      axis.text.x = if (show_x_axis_labels) element_text(size = 20, margin = margin(t = 6)) else element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 14, margin = margin(r = 4, l = 4)),
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      legend.position.inside = if (show_legend) c(0.98, 0.98) else NULL,
      legend.position = if (show_legend) "inside" else "none",
      legend.justification = if (show_legend) c(1, 1) else NULL,
      legend.background = if (show_legend) element_rect(fill = alpha("white", 0.7), color = NA) else NULL,
      legend.title = element_blank(),
      legend.text = if (show_legend) element_text(size = 17) else NULL,
      plot.margin = margin(t = if (reference_line == 0) 0 else 4, r = 6, b = if (reference_line == 0) 4 else 0, l = 4)
    )
  
  if (has_title) {
    # Dynamically push the title to the right ONLY if the panel has a y-axis title pushing the left boundary out
    p <- p + ggtitle(panel_title) + theme(
      plot.title = element_text(size = 20, hjust = 0, margin = margin(b = 10, l = if (show_y_axis_title) 45 else 0))
    )
  }
  
  return(p)
}