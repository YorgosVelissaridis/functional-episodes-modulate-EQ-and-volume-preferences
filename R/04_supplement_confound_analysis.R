# ==============================================================================
# 04_supplement_confounds.R: Confound Analysis and Volume-EQ Coupling
# ==============================================================================
source("R/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(brms)
  library(loo)
  library(posterior)
})

# ------------------------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------------------------
fit_or_load_model <- function(model_file, formula, family, data, force_refit) {
  if (file.exists(model_file) && !force_refit) {
    return(readRDS(model_file))
  }
  fit <- brm(
    formula = formula,
    data = data,
    family = family,
    backend = "cmdstanr",
    chains = CHAINS,
    iter = ITER,
    warmup = WARMUP,
    control = list(adapt_delta = ADAPT_DELTA),
    save_pars = save_pars(all = TRUE),
    file = NULL
  )
  saveRDS(fit, model_file)
  fit
}

loo_or_load <- function(loo_file, fit, force_recompute) {
  if (file.exists(loo_file) && !force_recompute) {
    return(readRDS(loo_file))
  }
  old_maxsize <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = FUTURE_GLOBALS_MAXSIZE)
  on.exit(options(future.globals.maxSize = old_maxsize), add = TRUE)
  
  loo_obj <- brms::loo(fit, moment_match = USE_MOMENT_MATCH, reloo = FALSE, cores = LOO_CORES)
  saveRDS(loo_obj, loo_file)
  loo_obj
}

compute_delta_from_loo <- function(loo_full, loo_target) {
  pw_full <- loo_full$pointwise[, "elpd_loo"]
  pw_target <- loo_target$pointwise[, "elpd_loo"]
  
  if (length(pw_full) != length(pw_target)) {
    stop("full and target models have different numbers of observations.")
  }
  
  d <- pw_target - pw_full
  N_obs <- length(d)
  delta_sum <- sum(d)
  se <- sqrt(N_obs * stats::var(d))
  
  tibble(
    delta_elpd = delta_sum,
    se = se,
    prob_multiplier = exp(delta_sum / N_obs)
  )
}

# ------------------------------------------------------------------------------
# 2. Data Preparation
# ------------------------------------------------------------------------------
cat("Preparing data for confound models...\n")
df <- read.csv(DATA_FILE) %>%
  mutate(
    excerpt_id = segment_id,
    selected_volume_db = selected_level_db,
    selected_eq = relevel(factor(selected_eq), ref = "Flat"),
    episode = relevel(factor(episode), ref = "Null"),
    noise = relevel(factor(noise), ref = "EV"),
    hearing_impairment = factor(hearing_impairment)
  )

# Standardize continuous variables for optimal MCMC sampling
df <- df %>%
  mutate(
    pre_arousal_z = as.numeric(scale(pre_arousal)),
    pre_valence_z = as.numeric(scale(pre_valence)),
    post_arousal_z = as.numeric(scale(post_arousal)),
    post_valence_z = as.numeric(scale(post_valence)),
    immersion_level_z = as.numeric(scale(immersion_level)),
    task_difficulty_z = as.numeric(scale(task_difficulty)),
    task_duration_sec_z = as.numeric(scale(task_duration_sec))
  )

# ------------------------------------------------------------------------------
# 3. Model Definition and Execution
# ------------------------------------------------------------------------------
confound_specs <- list(
  list(tag = "pre_mood", covars = "pre_arousal_z + pre_valence_z"),
  list(tag = "pre_post_mood", covars = "pre_arousal_z + pre_valence_z + post_arousal_z + post_valence_z"),
  list(tag = "hearing", covars = "hearing_impairment"),
  list(tag = "immersion", covars = "immersion_level_z"),
  list(tag = "difficulty", covars = "task_difficulty_z"),
  list(tag = "duration", covars = "task_duration_sec_z")
)

if (RUN_MODELING) {
  cat("\nLoading reference baseline models...\n")
  fit_eq <- readRDS(EQ_MODEL_FILE)
  fit_volume <- readRDS(VOLUME_MODEL_FILE)
  
  loo_eq <- loo_or_load(file.path(LOO_CACHE_DIR, "loo_full_eq.rds"), fit_eq, FALSE)
  loo_volume <- loo_or_load(file.path(LOO_CACHE_DIR, "loo_full_volume.rds"), fit_volume, FALSE)
  
  results_list <- list()
  
  # A. Test the 7 standard confounds
  for (spec in confound_specs) {
    cat(sprintf("\nEvaluating Confound: %s...\n", spec$tag))
    
    # EQ Model
    cat("  -> EQ Choice\n")
    eq_form <- update(EQ_FORMULA, paste(". ~ . +", spec$covars))
    fit_eq_conf <- fit_or_load_model(
      file.path(MODEL_CACHE_DIR, paste0("confound_eq_", spec$tag, ".rds")),
      eq_form, categorical(link = "logit"), df, FORCE_REFIT_CONFOUND_MODELS
    )
    loo_eq_conf <- loo_or_load(
      file.path(LOO_CACHE_DIR, paste0("loo_confound_eq_", spec$tag, ".rds")),
      fit_eq_conf, FORCE_RECOMPUTE_CONFOUND_LOO
    )
    res_eq <- compute_delta_from_loo(loo_eq, loo_eq_conf)
    res_eq$model <- "EQ"
    res_eq$confound <- spec$tag
    results_list[[length(results_list) + 1]] <- res_eq
    
    # Volume Model
    cat("  -> Volume Level\n")
    vol_form <- update(VOLUME_FORMULA, paste(". ~ . +", spec$covars))
    fit_vol_conf <- fit_or_load_model(
      file.path(MODEL_CACHE_DIR, paste0("confound_volume_", spec$tag, ".rds")),
      vol_form, gaussian(), df, FORCE_REFIT_CONFOUND_MODELS
    )
    loo_vol_conf <- loo_or_load(
      file.path(LOO_CACHE_DIR, paste0("loo_confound_volume_", spec$tag, ".rds")),
      fit_vol_conf, FORCE_RECOMPUTE_CONFOUND_LOO
    )
    res_vol <- compute_delta_from_loo(loo_volume, loo_vol_conf)
    res_vol$model <- "Volume"
    res_vol$confound <- spec$tag
    results_list[[length(results_list) + 1]] <- res_vol
  }
  
  # B. Test the specific Volume-predicting-EQ Model
  cat("\nEvaluating Confound: volume_predicting_eq...\n")
  cat("  -> EQ Choice\n")
  eq_vol_form <- update(EQ_FORMULA, . ~ . + selected_volume_db)
  fit_eq_vol <- fit_or_load_model(
    file.path(MODEL_CACHE_DIR, "confound_eq_volume.rds"),
    eq_vol_form, categorical(link = "logit"), df, FORCE_REFIT_CONFOUND_MODELS
  )
  loo_eq_vol <- loo_or_load(
    file.path(LOO_CACHE_DIR, "loo_confound_eq_volume.rds"),
    fit_eq_vol, FORCE_RECOMPUTE_CONFOUND_LOO
  )
  res_eq_vol <- compute_delta_from_loo(loo_eq, loo_eq_vol)
  res_eq_vol$model <- "EQ"
  res_eq_vol$confound <- "concurrent_volume"
  results_list[[length(results_list) + 1]] <- res_eq_vol
  
  # Save Master Confound Table
  confound_df <- bind_rows(results_list) %>%
    select(model, confound, delta_elpd, se, prob_multiplier) %>%
    arrange(model, confound)
  
  write.csv(confound_df, file.path(RESULTS_DIR, "confound_elpd_comparisons.csv"), row.names = FALSE)
  cat("\nConfound ELPD summary saved successfully.\n")
}

# ------------------------------------------------------------------------------
# 4. Volume-EQ Coupling Visualization
# ------------------------------------------------------------------------------
if (RUN_VISUALIZATION) {
  cat("\nGenerating Volume-EQ Coupling Plot...\n")
  
  EQ_COLOR_MAP <- c(
    "Flat" = "#7f7f7f", "BassCut" = "#fc8d62", "HighCut" = "#8da0cb",
    "VocalReduce" = "#e78ac3", "V-Shape" = "#a6d854", "VShape" = "#a6d854"
  )
  
  # Load the trained model
  fit_eq_vol <- readRDS(file.path(MODEL_CACHE_DIR, "confound_eq_volume.rds"))
  
  # Save the sample size diagnostics
  under_minus10_counts <- df %>%
    filter(!is.na(selected_volume_db), selected_volume_db <= -10) %>%
    count(selected_volume_db, name = "sample_size") %>%
    arrange(selected_volume_db)
  
  under_minus10_total <- under_minus10_counts %>%
    summarise(sample_size = sum(sample_size, na.rm = TRUE)) %>%
    mutate(selected_volume_db = "Total (<= -10 dB)") %>%
    select(selected_volume_db, sample_size)
  
  under_minus10_out <- bind_rows(
    under_minus10_counts %>% mutate(selected_volume_db = as.character(selected_volume_db)),
    under_minus10_total
  )
  write.csv(under_minus10_out, file.path(RESULTS_DIR, "full_plus_volume_eq_under_minus10_sample_sizes.csv"), row.names = FALSE)

  # Compute posterior odds ratios
  volume_grid <- seq(floor(min(df$selected_volume_db, na.rm = TRUE)), ceiling(max(df$selected_volume_db, na.rm = TRUE)), by = 1)
  newdata_grid <- data.frame(
    selected_volume_db = volume_grid,
    noise = factor("EV", levels = levels(df$noise)),
    episode = factor("Null", levels = levels(df$episode))
  )
  newdata_baseline <- data.frame(
    selected_volume_db = 0,
    noise = factor("EV", levels = levels(df$noise)),
    episode = factor("Null", levels = levels(df$episode))
  )
  
  epred_grid <- posterior_epred(fit_eq_vol, newdata = newdata_grid, re_formula = NA)
  epred_0db <- posterior_epred(fit_eq_vol, newdata = newdata_baseline, re_formula = NA)[, 1, ]
  
  eq_levels <- dimnames(epred_grid)[[3]]
  if (is.null(eq_levels)) eq_levels <- levels(df$selected_eq)
  flat_idx <- match("Flat", eq_levels)
  
  eps <- 1e-8
  plot_eq_levels <- setdiff(eq_levels, "Flat")
  
  plot_df <- bind_rows(lapply(plot_eq_levels, function(eq_name) {
    k <- match(eq_name, eq_levels)
    p_k_v <- pmax(epred_grid[, , k], eps)
    p_flat_v <- pmax(epred_grid[, , flat_idx], eps)
    p_k_0 <- pmax(epred_0db[, k], eps)
    p_flat_0 <- pmax(epred_0db[, flat_idx], eps)
    
    odds_vs_flat_v <- p_k_v / p_flat_v
    odds_vs_flat_0 <- p_k_0 / p_flat_0
    or_mat <- odds_vs_flat_v / odds_vs_flat_0
    
    data.frame(
      selected_volume_db = volume_grid,
      selected_eq = eq_name,
      odds_ratio = apply(or_mat, 2, mean),
      lower = apply(or_mat, 2, quantile, probs = 0.025),
      upper = apply(or_mat, 2, quantile, probs = 0.975),
      stringsAsFactors = FALSE
    )
  }))
  
  write.csv(plot_df, file.path(RESULTS_DIR, "full_plus_volume_eq.csv"), row.names = FALSE)
  
  # Generate Plot
  x_limits <- range(df$selected_volume_db, na.rm = TRUE)
  hist_counts <- hist(df$selected_volume_db, breaks = seq(floor(x_limits[1]), ceiling(x_limits[2]) + 1, by = 1), plot = FALSE)$counts
  hist_max <- max(hist_counts, na.rm = TRUE)
  if (!is.finite(hist_max) || hist_max <= 0) hist_max <- 1
  
  y_min <- 0
  y_max <- 20
  x_breaks <- seq(floor(x_limits[1] / 5) * 5, ceiling(x_limits[2] / 5) * 5, by = 5)
  y_breaks <- c(0, 1, 3, 6, 9, 12, 15, 18)
  
  p_combined <- ggplot() +
    geom_histogram(
      data = df,
      aes(x = selected_volume_db, y = after_stat((count / hist_max) * y_max)),
      binwidth = 1, center = 0, closed = "left",
      fill = "#5a5a5a", color = "white", linewidth = 0.25, alpha = 0.52
    ) +
    geom_ribbon(
      data = plot_df,
      aes(x = selected_volume_db, ymin = lower, ymax = upper, fill = selected_eq),
      alpha = 0.16, linewidth = 0, show.legend = FALSE
    ) +
    geom_line(
      data = plot_df,
      aes(x = selected_volume_db, y = odds_ratio, color = selected_eq),
      linewidth = 1
    ) +
    geom_point(
      data = plot_df,
      aes(x = selected_volume_db, y = odds_ratio, color = selected_eq),
      size = 1.6
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8, color = "grey40") +
    coord_cartesian(xlim = x_limits, ylim = c(y_min, y_max)) +
    scale_x_continuous(breaks = x_breaks) +
    scale_color_manual(values = EQ_COLOR_MAP, drop = TRUE) +
    scale_fill_manual(values = EQ_COLOR_MAP, drop = TRUE) +
    scale_y_continuous(
      name = "Odds Ratio (vs Flat, baseline 0 dB)",
      breaks = y_breaks,
      sec.axis = sec_axis(~ (. / y_max) * hist_max, name = "Number of responses")
    ) +
    labs(x = "Volume Level [dB]", color = "EQ") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
      legend.position = "inside",
      legend.position.inside = c(0.18, 0.97),
      legend.justification = c(0, 1)
    )
  
  out_fig_path <- file.path(RESULTS_DIR, "full_plus_volume_eq.pdf")
  ggsave(out_fig_path, p_combined, width = 8.5, height = 4.8)
  cat("Plot saved successfully to:", out_fig_path, "\n")
}