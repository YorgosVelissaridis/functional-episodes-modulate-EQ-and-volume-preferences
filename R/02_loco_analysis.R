# ==============================================================================
# 02_loco_analysis.R: Leave-One-Covariate-Out (LOCO) ELPD Analysis
# ==============================================================================

source("R/00_config.R")

# ------------------------------
# Data Preparation
# ------------------------------
prepare_model_data <- function(path) {
  df <- read.csv(path)
  
  # Make octave-band column names formula-safe.
  if ("oct_level_db_31.5" %in% names(df) && !"oct_level_db_31_5" %in% names(df)) {
    df <- dplyr::rename(df, oct_level_db_31_5 = `oct_level_db_31.5`)
  }
  
  df %>%
    dplyr::rename(excerpt_id = segment_id, selected_volume_db = selected_level_db) %>%
    mutate(
      gender_merged = ifelse(
        tolower(trimws(as.character(gender))) %in% c("man", "women"),
        tolower(trimws(as.character(gender))),
        "others"
      ) %>% factor(levels = c("man", "women", "others")),
      genres = if ("genres" %in% names(df)) factor(genres) else NULL,
      selected_eq = relevel(factor(selected_eq), ref = "Flat"),
      episode = relevel(factor(episode), ref = "Null"),
      noise = relevel(factor(noise), ref = "EV")
    )
}

# ------------------------------
# Modeling and LOO Helpers
# ------------------------------
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
  
  compute_loo <- function(use_reloo) {
    old_maxsize <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = FUTURE_GLOBALS_MAXSIZE)
    on.exit(options(future.globals.maxSize = old_maxsize), add = TRUE)
    
    suppressWarnings(
      brms::loo(
        fit,
        moment_match = USE_MOMENT_MATCH,
        reloo = use_reloo,
        cores = LOO_CORES
      )
    )
  }
  
  loo_obj <- compute_loo(RELOO)
  if (!RELOO && AUTO_RELOO_IF_HIGH_K) {
    has_high_k <- tryCatch({
      any(loo::pareto_k_values(loo_obj) > 0.7, na.rm = TRUE)
    }, error = function(e) FALSE)
    
    if (has_high_k) {
      cat("  [High Pareto k detected. Refitting problematic points exactly via reloo...]\n")
      loo_obj <- compute_loo(TRUE)
    }
  }
  
  saveRDS(loo_obj, loo_file)
  loo_obj
}

#' Compute out-of-sample predictive difference (Delta ELPD) between two models.
#'
#' Calculates the exact pointwise ELPD drop when a covariate is removed. 
#' Includes a Bayesian Bootstrap to estimate the standard error and 95% 
#' credible intervals, capturing the asymmetric skew of the predictive distributions.
#' 
#' @param loo_full The loo object for the full baseline model.
#' @param loo_target The loo object for the reduced model.
#' @return A tibble containing the ELPD difference, SE, and bootstrap quantiles.
compute_delta_from_loo <- function(loo_full, loo_target) {
  pw_full <- loo_full$pointwise[, "elpd_loo"]
  pw_target <- loo_target$pointwise[, "elpd_loo"]
  
  if (length(pw_full) != length(pw_target)) {
    stop("full and target models have different numbers of observations.")
  }
  
  d <- pw_target - pw_full
  N_obs <- length(d)
  delta_sum <- sum(d)
  
  # Standard Error based on variance of pointwise differences
  se <- sqrt(N_obs * stats::var(d))
  
  # Bayesian Bootstrap for the summed ELPD difference
  B <- 4000
  bb_sums <- numeric(B)
  for (i in seq_len(B)) {
    w <- rexp(N_obs, 1)        # Draw standard exponentials
    w <- w / mean(w)           # Normalize so weights sum to N_obs
    bb_sums[i] <- sum(w * d)   # Bootstrapped sum
  }
  
  # Extract the true asymmetric 2.5% and 97.5% quantiles
  q_2_5 <- quantile(bb_sums, 0.025)
  q_97_5 <- quantile(bb_sums, 0.975)
  
  tibble(
    delta_elpd = delta_sum,
    se = se,
    N = N_obs,
    delta_elpd_per_obs = delta_sum / N_obs,
    prob_multiplier = exp(delta_sum / N_obs),
    ymin = q_2_5,
    ymax = q_97_5
  )
}

# ------------------------------
# Plot Data Helpers
# ------------------------------
build_plot_df <- function(delta_df) {
  plot_rows <- list()
  idx <- 1L
  
  for (metric_name in c("EQ", "Volume")) {
    for (grp in PLOT_GROUPS) {
      grp_df <- delta_df %>%
        filter(metric == metric_name, model %in% grp$model) %>%
        mutate(
          panel = grp$panel,
          model = factor(model, levels = grp$model)
        ) %>%
        arrange(model)
      
      plot_rows[[idx]] <- grp_df
      idx <- idx + 1L
    }
  }
  
  bind_rows(plot_rows) %>%
    mutate(
      metric = factor(metric, levels = c("EQ", "Volume")),
      panel = factor(panel, levels = vapply(PLOT_GROUPS, function(x) x$panel, character(1)))
    )
}

save_delta_plot <- function(plot_df, out_pdf) {
  panel_levels <- vapply(PLOT_GROUPS, function(x) x$panel, character(1))
  
  make_panel_plot <- function(
    df_panel,
    panel_title,
    y_limits = NULL,
    y_label = NULL,
    show_xticks = TRUE,
    show_y_axis = TRUE,
    margin_top = 6,
    margin_bottom = 2
  ) {
    x_levels <- levels(df_panel$model)
    x_labels <- x_levels
    names(x_labels) <- x_levels
    has_manual <- x_levels %in% names(MODEL_X_LABELS)
    x_labels[has_manual] <- unname(MODEL_X_LABELS[x_levels[has_manual]])
    
    p <- ggplot(df_panel, aes(x = model, y = delta_elpd)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, linewidth = 0.6) +
      geom_point(size = 2.6) +
      scale_x_discrete(labels = x_labels) +
      labs(title = panel_title, x = NULL, y = y_label) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(
          angle = 0, hjust = 0.5, vjust = 0.5,
          size = AXIS_TICK_TEXT_SIZE, color = AXIS_TICK_TEXT_COLOR
        ),
        axis.text.y = element_text(size = AXIS_TICK_TEXT_SIZE, color = AXIS_TICK_TEXT_COLOR),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = PLOT_TITLE_SIZE, face = "plain", hjust = 0.5),
        plot.margin = margin(t = margin_top, r = 2, b = margin_bottom, l = 2)
      )
    
    if (!show_xticks) {
      p <- p + theme(
        axis.text.x = element_text(color = "transparent", angle = 0, hjust = 0.5, vjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(t = margin_top, r = 2, b = margin_bottom, l = 2)
      )
    }
    
    if (!show_y_axis) {
      p <- p + theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()
      )
    }
    
    if (!is.null(y_limits)) {
      p <- p + coord_cartesian(ylim = y_limits)
    }
    
    p
  }
  
  eq_bridge_y <- plot_df %>%
    filter(metric == "EQ", panel == panel_levels[1], model == "drop_excerpt") %>%
    summarise(y = mean(delta_elpd, na.rm = TRUE)) %>%
    pull(y)
  if (!is.finite(eq_bridge_y)) eq_bridge_y <- NA_real_
  
  eq_bridge_y_participant <- plot_df %>%
    filter(metric == "EQ", panel == panel_levels[1], model == "drop_participant") %>%
    summarise(y = mean(delta_elpd, na.rm = TRUE)) %>%
    pull(y)
  if (!is.finite(eq_bridge_y_participant)) eq_bridge_y_participant <- NA_real_
  
  volume_bridge_y_participant <- plot_df %>%
    filter(metric == "Volume", panel == panel_levels[1], model == "drop_participant") %>%
    summarise(y = mean(delta_elpd, na.rm = TRUE)) %>%
    pull(y)
  if (!is.finite(volume_bridge_y_participant)) volume_bridge_y_participant <- NA_real_
  
  volume_min <- min(plot_df$ymin[plot_df$metric == "Volume"], na.rm = TRUE)
  volume_max <- max(plot_df$ymax[plot_df$metric == "Volume"], na.rm = TRUE)
  volume_limits <- c(volume_min, volume_max)
  
  eq_plots <- lapply(seq_along(panel_levels), function(i) {
    pnl <- panel_levels[i]
    df_panel <- plot_df %>% filter(metric == "EQ", panel == pnl)
    
    make_panel_plot(
      df_panel,
      pnl,
      y_limits = c(-170, 0),
      y_label = if (i == 1) Y_AXIS_LABEL_EXPR else NULL,
      show_xticks = TRUE,
      show_y_axis = (i == 1)
    )
  })
  
  volume_a <- make_panel_plot(
    plot_df %>% filter(metric == "Volume", panel == panel_levels[1]),
    "(d) Impact of Factor Removal\non Volume Prediction",
    y_limits = volume_limits,
    y_label = Y_AXIS_LABEL_EXPR,
    show_y_axis = TRUE,
    margin_top = 6,
    margin_bottom = 2
  )
  volume_blank <- ggplot() + theme_void()
  volume_c <- make_panel_plot(
    plot_df %>% filter(metric == "Volume", panel == panel_levels[3]),
    "(e) Predictive Contribution of Participant Traits\non Volume Prediction",
    y_limits = volume_limits,
    y_label = NULL,
    show_y_axis = FALSE,
    margin_top = 6,
    margin_bottom = 2
  )
  volume_plots <- list(volume_a, volume_blank, volume_c)
  
  draw_layout <- function() {
    eq_grobs <- lapply(eq_plots, ggplotGrob)
    volume_grobs <- lapply(volume_plots, ggplotGrob)
    
    for (i in seq_along(eq_grobs)) {
      max_heights <- grid::unit.pmax(eq_grobs[[i]]$heights, volume_grobs[[i]]$heights)
      eq_grobs[[i]]$heights <- max_heights
      volume_grobs[[i]]$heights <- max_heights
    }
    
    lay <- grid::grid.layout(
      nrow = 3,
      ncol = 5,
      widths = grid::unit.c(
        grid::unit(1, "null"),
        grid::unit(SUBPLOT_COL_SPACER_IN, "in"),
        grid::unit(1, "null"),
        grid::unit(SUBPLOT_COL_SPACER_IN, "in"),
        grid::unit(1, "null")
      ),
      heights = grid::unit.c(
        grid::unit(1, "null"),
        grid::unit(0.28, "in"),
        grid::unit(1, "null")
      )
    )
    
    grid::pushViewport(grid::viewport(layout = lay))
    panel_cols <- c(1, 3, 5)
    
    for (i in 1:3) {
      grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = panel_cols[i]))
      grid::grid.draw(eq_grobs[[i]])
      grid::upViewport()
    }
    
    for (i in 1:3) {
      grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = panel_cols[i]))
      grid::grid.draw(volume_grobs[[i]])
      grid::upViewport()
    }
    
    # Dashed links are drawn globally so they are not clipped by panel boundaries.
    if (is.finite(eq_bridge_y) || is.finite(eq_bridge_y_participant)) {
      grid::pushViewport(
        grid::viewport(
          layout.pos.row = 1,
          layout.pos.col = 1:5,
          xscale = c(1, 12),
          yscale = c(-170, 0)
        )
      )
      if (is.finite(eq_bridge_y)) {
        grid::grid.lines(
          x = grid::unit(c(3.52, 5.24), "native"),
          y = grid::unit(c(-87, -87), "native"),
          gp = grid::gpar(col = "grey65", lwd = 3.0, lty = "dashed")
        )
      }
      if (is.finite(eq_bridge_y_participant)) {
        grid::grid.lines(
          x = grid::unit(c(4.22, 8.92), "native"),
          y = grid::unit(c(-109, -109), "native"),
          gp = grid::gpar(col = "grey65", lwd = 3.0, lty = "dashed")
        )
      }
      grid::upViewport()
    }
    
    if (is.finite(volume_bridge_y_participant)) {
      grid::pushViewport(
        grid::viewport(
          layout.pos.row = 3,
          layout.pos.col = 1:5,
          xscale = c(1, 12),
          yscale = volume_limits
        )
      )
      grid::grid.lines(
        x = grid::unit(c(4.2, 8.9), "native"),
        y = grid::unit(c(-650, -650), "native"),
        gp = grid::gpar(col = "grey65", lwd = 3.0, lty = "dashed")
      )
      grid::upViewport()
    }
    
    grid::upViewport()
  }
  
  grDevices::pdf(OUT_ELPD_PLOT, width = 13, height = 9)
  draw_layout()
  grDevices::dev.off()
}

# ------------------------------
# Cache Utilities
# ------------------------------
all_model_tags <- vapply(MODEL_SPECS, function(x) x$tag, character(1))

has_delta_row <- function(df_in, metric_name, tag) {
  any(df_in$metric == metric_name & df_in$model == tag)
}

# ------------------------------
# Main Execution
# ------------------------------
df <- prepare_model_data(DATA_FILE)

delta_df <- NULL
if (file.exists(DELTA_CACHE_CSV) && !FORCE_RECOMPUTE_DELTA) {
  delta_df <- read.csv(DELTA_CACHE_CSV, stringsAsFactors = FALSE)
}
if (is.null(delta_df)) {
  delta_df <- tibble(
    metric = character(), 
    model = character(), 
    delta_elpd = numeric(), 
    se = numeric(),
    N = numeric(), 
    delta_elpd_per_obs = numeric(), 
    prob_multiplier = numeric(),
    ymin = numeric(),
    ymax = numeric()
  )
}

if (RUN_MODELING) {
  full_refit <- FORCE_REFIT_FULL_MODELS
  full_loo_recompute <- FORCE_RECOMPUTE_FULL_LOO || full_refit
  
  full_eq <- fit_or_load_model(
    EQ_MODEL_FILE,
    EQ_FORMULA,
    categorical(link = "logit"),
    df,
    full_refit
  )
  full_volume <- fit_or_load_model(
    VOLUME_MODEL_FILE,
    VOLUME_FORMULA,
    gaussian(),
    df,
    full_refit
  )
  
  loo_full_eq <- loo_or_load(file.path(LOO_CACHE_DIR, "loo_full_eq.rds"), full_eq, full_loo_recompute)
  loo_full_volume <- loo_or_load(file.path(LOO_CACHE_DIR, "loo_full_volume.rds"), full_volume, full_loo_recompute)
  
  for (spec in MODEL_SPECS) {
    tag <- spec$tag
    force_refit_tag <- FORCE_REFIT_LOCO_MODELS
    force_loo_tag <- FORCE_RECOMPUTE_LOCO_LOO || force_refit_tag
    force_delta_tag <- FORCE_RECOMPUTE_DELTA
    
    missing_cache <- !has_delta_row(delta_df, "EQ", tag) || !has_delta_row(delta_df, "Volume", tag)
    needs_update <- force_delta_tag || force_loo_tag || missing_cache
    if (!needs_update) {
      next
    }
    
    target_eq <- fit_or_load_model(
      file.path(MODEL_CACHE_DIR, paste0(tag, "_eq.rds")),
      spec$eq,
      categorical(link = "logit"),
      df,
      force_refit_tag
    )
    target_volume <- fit_or_load_model(
      file.path(MODEL_CACHE_DIR, paste0(tag, "_volume.rds")),
      spec$volume,
      gaussian(),
      df,
      force_refit_tag
    )
    
    loo_target_eq <- loo_or_load(
      file.path(LOO_CACHE_DIR, paste0("loo_", tag, "_eq.rds")),
      target_eq,
      force_loo_tag
    )
    loo_target_volume <- loo_or_load(
      file.path(LOO_CACHE_DIR, paste0("loo_", tag, "_volume.rds")),
      target_volume,
      force_loo_tag
    )
    
    eq_delta <- compute_delta_from_loo(loo_full_eq, loo_target_eq)
    volume_delta <- compute_delta_from_loo(loo_full_volume, loo_target_volume)
    
    # Store the new columns for both metrics
    delta_df <- delta_df %>%
      filter(!(model == tag & metric %in% c("EQ", "Volume"))) %>%
      bind_rows(
        tibble(
          metric = "EQ", model = tag, 
          delta_elpd = eq_delta$delta_elpd, se = eq_delta$se,
          N = eq_delta$N, delta_elpd_per_obs = eq_delta$delta_elpd_per_obs, prob_multiplier = eq_delta$prob_multiplier,
          ymin = eq_delta$ymin, ymax = eq_delta$ymax
        ),
        tibble(
          metric = "Volume", model = tag, 
          delta_elpd = volume_delta$delta_elpd, se = volume_delta$se,
          N = volume_delta$N, delta_elpd_per_obs = volume_delta$delta_elpd_per_obs, prob_multiplier = volume_delta$prob_multiplier,
          ymin = volume_delta$ymin, ymax = volume_delta$ymax
        )
      )
  }
  
  delta_df <- delta_df %>%
    mutate(
      metric = factor(metric, levels = c("EQ", "Volume")),
      model = factor(model, levels = all_model_tags)
    ) %>%
    arrange(metric, model) %>%
    mutate(
      metric = as.character(metric),
      model = as.character(model)
    )
  
  write.csv(delta_df, DELTA_CACHE_CSV, row.names = FALSE)
  saveRDS(delta_df, DELTA_CACHE_RDS)
}

if (nrow(delta_df) == 0) {
  if (file.exists(DELTA_CACHE_RDS)) {
    delta_df <- readRDS(DELTA_CACHE_RDS)
  } else {
    stop("Delta ELPD cache not found. Run with RUN_MODELING = TRUE first.")
  }
}

if (RUN_LOCO_PLOTTING) {
  plot_df <- build_plot_df(delta_df)
  save_delta_plot(plot_df, OUT_ELPD_PLOT)
}

# -------------------
# Table Export
# -------------------
cat("\nExporting Interpretable ELPD Reductions table...\n")

elpd_table <- delta_df %>%
  mutate(
    delta_elpd = round(delta_elpd, 2),
    se = round(se, 2),
    ymin = round(ymin, 2),
    ymax = round(ymax, 2),
    delta_elpd_per_obs = round(delta_elpd_per_obs, 4),
    prob_multiplier = round(prob_multiplier, 4),
    pct_drop = paste0(round((1 - prob_multiplier) * 100, 1), "%")
  ) %>%
  select(metric, model, delta_elpd, ymin, ymax, se, delta_elpd_per_obs, prob_multiplier, pct_drop)

write.csv(elpd_table, OUT_ELPD_TABLE, row.names = FALSE)
cat("ELPD reductions table exported successfully to: ", OUT_ELPD_TABLE, "\n")