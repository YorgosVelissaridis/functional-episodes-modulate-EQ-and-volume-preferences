# ==============================================================================
# 03_supplement_sensitivity.R: Individual Episode Sensitivity Analysis
# Calculates within-subject Volume SD and EQ Perplexity across functional episodes
# ==============================================================================
source("R/00_config.R")

# 1. Setup and Libraries -------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
})

# 2. Mathematical Helper Functions ---------------------------------------------
# Calculate Shannon Entropy in bits (Base 2)
calc_entropy <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  
  props <- prop.table(table(x))
  props <- props[props > 0] 
  
  -sum(props * log2(props))
}

# 3. Main Data Processing ------------------------------------------------------
cat("Loading and processing data for sensitivity analysis...\n")
df <- read.csv(DATA_FILE, stringsAsFactors = FALSE) %>%
  dplyr::rename(selected_volume_db = selected_level_db)

# Step A: Calculate metrics per participant, per noise condition
df_condition_sens <- df %>%
  group_by(participant_id, noise) %>%
  summarise(
    vol_sd = sd(selected_volume_db, na.rm = TRUE),
    eq_entropy = calc_entropy(selected_eq),
    .groups = "drop"
  )

# Step B: Average the metrics across EV and Diesel, then convert to Perplexity
df_participant_sens <- df_condition_sens %>%
  group_by(participant_id) %>%
  summarise(
    mean_vol_sd = mean(vol_sd, na.rm = TRUE),
    mean_eq_entropy = mean(eq_entropy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    eq_perplexity = 2^mean_eq_entropy
  ) %>%
  filter(!is.na(mean_vol_sd) & !is.na(eq_perplexity))

# 4. Generate Descriptive Statistics -------------------------------------------
cat("\n--- Descriptive Statistics ---\n")
cat("EQ Perplexity:\n")
cat("Median:", median(df_participant_sens$eq_perplexity), "\n")
cat("Min:   ", min(df_participant_sens$eq_perplexity), "\n")
cat("Max:   ", max(df_participant_sens$eq_perplexity), "\n")

purists <- df_participant_sens %>% filter(eq_perplexity == 1.0)
cat("Number of EQ Purists (Perplexity = 1.0):", nrow(purists), "\n")

cat("\nVolume SD:\n")
cat("Median:", median(df_participant_sens$mean_vol_sd), "dB\n")
cat("Min:   ", min(df_participant_sens$mean_vol_sd), "dB\n")
cat("Max:   ", max(df_participant_sens$mean_vol_sd), "dB\n")

reactive_listeners <- df_participant_sens %>% filter(mean_vol_sd >= 5.0)
cat("Number of highly reactive volume listeners (SD >= 5.0):", nrow(reactive_listeners), "\n")

# 5. Plotting ------------------------------------------------------------------
if (RUN_VISUALIZATION) {
  cat("\nGenerating histograms...\n")
  
  plot_a <- ggplot(df_participant_sens, aes(x = eq_perplexity)) +
    geom_histogram(binwidth = 0.05, fill = "#4c8c5c", color = "black", alpha = 0.8) +
    labs(
      title = "(a) Perplexity of Preferred EQ Profiles",
      subtitle = "Average of EV & Diesel Trials",
      x = "Perplexity (Effective EQ Count)",
      y = "Number of Participants"
    ) +
    scale_x_continuous(breaks = seq(1, 5, by = 1), limits = c(0.8, 5.2)) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0, face = "plain"))
  
  plot_b <- ggplot(df_participant_sens, aes(x = mean_vol_sd)) +
    geom_histogram(binwidth = 0.5, fill = "#5b9bd5", color = "black", alpha = 0.8) +
    labs(
      title = "(b) SD of Preferred Volume Level",
      subtitle = "Average of EV & Diesel Trials",
      x = "Volume SD (dB)",
      y = ""
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0, face = "plain"))
  
  combined_plot <- gridExtra::grid.arrange(plot_a, plot_b, ncol = 2)
  
  out_plot_path <- file.path(RESULTS_DIR, "episode_sensitivity_histograms.pdf")
  ggsave(out_plot_path, combined_plot, width = 12, height = 5)
  cat("Plot saved successfully to:", out_plot_path, "\n")
}