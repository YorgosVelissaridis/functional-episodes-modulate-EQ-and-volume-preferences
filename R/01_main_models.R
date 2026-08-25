# ==============================================================================
# 01_main_models.R: Main Hierarchical Models - Fitting, R2 Extraction, Tables, and Visualization
# ==============================================================================

source("R/00_config.R")

# 1. Load and Encode Data ------------------------------------------------------
model_df <- read.csv(DATA_FILE) %>%
  rename(excerpt_id = segment_id, selected_volume_db = selected_level_db) %>%
  mutate(
    selected_eq = relevel(factor(selected_eq), ref = "Flat"),
    episode = relevel(factor(episode), ref = "Null"),
    noise = relevel(factor(noise), ref = "EV")
  )

# 2. Fit or Load Models --------------------------------------------------------
if (RUN_MODELING && (!file.exists(EQ_MODEL_FILE) || FORCE_REFIT_EQ)) {
  fit_eq <- brm(
    formula = EQ_FORMULA,
    data = model_df,
    family = categorical(link = "logit"),
    backend = "cmdstanr",
    chains = CHAINS,
    iter = ITER,
    warmup = WARMUP,
    control = list(adapt_delta = ADAPT_DELTA),
    save_pars = save_pars(all = TRUE),
    file = NULL
  )
  saveRDS(fit_eq, EQ_MODEL_FILE)
} else {
  fit_eq <- readRDS(EQ_MODEL_FILE)
}

if (RUN_MODELING && (!file.exists(VOLUME_MODEL_FILE) || FORCE_REFIT_VOLUME)) {
  fit_volume <- brm(
    formula = VOLUME_FORMULA,
    data = model_df,
    family = gaussian(),
    backend = "cmdstanr",
    chains = CHAINS,
    iter = ITER,
    warmup = WARMUP,
    control = list(adapt_delta = ADAPT_DELTA),
    save_pars = save_pars(all = TRUE),
    file = NULL
  )
  saveRDS(fit_volume, VOLUME_MODEL_FILE)
} else {
  fit_volume <- readRDS(VOLUME_MODEL_FILE)
}

# 3. LOO-Adjusted Pseudo-R2 Calculation ----------------------------------------
cat("\nComputing LOO-Adjusted Pseudo-R2 Metrics...\n")

r2_volume_file <- file.path(LOO_CACHE_DIR, "loo_r2_volume_main.rds")
if (file.exists(r2_volume_file) && !FORCE_RECOMPUTE_FULL_LOO) {
  r2_volume <- readRDS(r2_volume_file)
} else {
  r2_volume <- loo_R2(fit_volume, cores = LOO_CORES)
  saveRDS(r2_volume, r2_volume_file)
}

r2_volume_df <- as_tibble(r2_volume) %>% 
  mutate(Model = "Volume (Gaussian - LOO-R2)") %>%
  select(Model, Estimate, Est.Error, Q2.5, Q97.5)

if (RUN_MODELING && (!file.exists(NULL_EQ_MODEL_FILE) || FORCE_REFIT_EQ)) {
  fit_eq_null <- brm(
    formula = selected_eq ~ 1,
    data = model_df,
    family = categorical(link = "logit"),
    backend = "cmdstanr",
    chains = CHAINS,
    iter = ITER,
    warmup = WARMUP,
    file = NULL
  )
  saveRDS(fit_eq_null, NULL_EQ_MODEL_FILE)
} else {
  fit_eq_null <- readRDS(NULL_EQ_MODEL_FILE)
}

loo_eq_full_file <- file.path(LOO_CACHE_DIR, "loo_eq_full_main.rds")
if (file.exists(loo_eq_full_file) && !FORCE_RECOMPUTE_FULL_LOO) {
  loo_eq_full <- readRDS(loo_eq_full_file)
} else {
  loo_eq_full <- loo(fit_eq, cores = LOO_CORES)
  saveRDS(loo_eq_full, loo_eq_full_file)
}

loo_eq_null_file <- file.path(LOO_CACHE_DIR, "loo_eq_null_main.rds")
if (file.exists(loo_eq_null_file) && !FORCE_RECOMPUTE_FULL_LOO) {
  loo_eq_null <- readRDS(loo_eq_null_file)
} else {
  loo_eq_null <- loo(fit_eq_null, cores = LOO_CORES)
  saveRDS(loo_eq_null, loo_eq_null_file)
}

mcfadden_r2 <- 1 - (loo_eq_full$estimates["elpd_loo", "Estimate"] / loo_eq_null$estimates["elpd_loo", "Estimate"])

r2_eq_df <- tibble(
  Model = "EQ Profile (Categorical - LOO McFadden)",
  Estimate = mcfadden_r2,
  Est.Error = NA_real_, 
  Q2.5 = NA_real_,
  Q97.5 = NA_real_
)

r2_table <- bind_rows(r2_eq_df, r2_volume_df) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write.csv(r2_table, OUT_R2_TABLE, row.names = FALSE)
cat("Out-of-sample LOO-R2 summary saved successfully.\n")

# 4. Posterior Processing for Tables and Plots ---------------------------------
cat("\nExtracting and summarizing posteriors...\n")
source("R/utils.R")

# Posterior Processing: EQ Model
coef_eq <- as_draws_df(fit_eq)
fixed_terms_eq <- names(coef_eq)[grepl("^b_mu", names(coef_eq)) & grepl("_(episode|noise)", names(coef_eq))]

eq_fixed_draws <- map_dfr(fixed_terms_eq, function(term) {
  condition <- if (grepl("_episode", term)) sub("^.*_episode", "", term) else sub("^.*_noise", "", term)
  tibble(
    draw = coef_eq$.draw,
    condition = recode(condition, "FocusMMotivation" = "Focus-Motivation", .default = condition),
    eq = sub("^b_mu", "", sub("_(episode|noise).*", "", term)),
    or = exp(coef_eq[[term]])
  )
}) %>%
  filter(eq != "Flat")

sd_terms_eq <- names(coef_eq)[grepl("^sd_(participant_id|excerpt_id)__mu", names(coef_eq)) & grepl("_Intercept$", names(coef_eq))]
mor_k <- sqrt(2) * qnorm(0.75)

eq_random_draws <- map_dfr(sd_terms_eq, function(term) {
  group_name <- sub("^sd_(participant_id|excerpt_id)__.*$", "\\1", term)
  eq_name <- sub("_Intercept$", "", sub("^sd_(participant_id|excerpt_id)__mu", "", term))
  tibble(
    condition = ifelse(group_name == "excerpt_id", "Excerpt", "Individual"),
    eq = eq_name,
    or = exp(mor_k * coef_eq[[term]])
  )
}) %>%
  filter(eq != "Flat")

eq_fixed_summary <- summarise_or(eq_fixed_draws, "or", c("condition", "eq"))
eq_random_summary <- summarise_or(eq_random_draws, "or", c("condition", "eq"))

# Posterior Processing: Volume Model
coef_volume <- as_draws_df(fit_volume)
fixed_terms_volume <- names(coef_volume)[grepl("^b_", names(coef_volume)) & grepl("_(episode|noise)", names(coef_volume))]

volume_fixed_draws <- map_dfr(fixed_terms_volume, function(term) {
  condition <- if (grepl("_episode", term)) sub("^.*_episode", "", term) else sub("^.*_noise", "", term)
  tibble(
    condition = recode(condition, "FocusMMotivation" = "Focus-Motivation", .default = condition),
    eq = "Volume",
    or = coef_volume[[term]]
  )
})

volume_random_draws <- bind_rows(
  map_dfr(c("participant_id", "excerpt_id"), function(group_name) {
    term <- paste0("sd_", group_name, "__Intercept")
    tibble(
      condition = ifelse(group_name == "excerpt_id", "Excerpt", "Individual"),
      eq = "Volume",
      or = coef_volume[[term]]
    )
  }),
  tibble(
    condition = "Residual (Sigma)",
    eq = "Volume",
    or = coef_volume$sigma
  )
)

volume_fixed_summary <- summarise_or(volume_fixed_draws, "or", c("condition", "eq"))
volume_random_summary <- summarise_or(volume_random_draws, "or", c("condition", "eq"))

# 5. Table Formatting and Export -----------------------------------------------
cat("Exporting summary tables...\n")
factor_order <- c(
  factor_label_for_table(c("Enjoyment", "Distraction", "Relaxation", "Focus-Motivation")),
  factor_label_for_table(c("DIESEL")),
  factor_label_for_table(c("Excerpt", "Individual"))
)

eq_table <- bind_rows(
  summarise_plot_stats(eq_fixed_draws) %>%
    transmute(
      Factor = factor_label_for_table(condition),
      EQ_Profile = eq,
      Odds_Ratio = mean,
      Median_Odds_Ratio = NA_real_,
      SD = sd,
      CI_2_5 = q2_5,
      CI_97_5 = q97_5
    ),
  summarise_plot_stats(eq_random_draws) %>%
    transmute(
      Factor = factor_label_for_table(condition),
      EQ_Profile = eq,
      Odds_Ratio = NA_real_,
      Median_Odds_Ratio = median,
      SD = sd,
      CI_2_5 = q2_5,
      CI_97_5 = q97_5
    )
) %>%
  arrange(match(Factor, factor_order), EQ_Profile) %>%
  mutate(across(c(Odds_Ratio, Median_Odds_Ratio, SD, CI_2_5, CI_97_5), ~ round(.x, 2)))

volume_table <- bind_rows(
  summarise_plot_stats(volume_fixed_draws) %>%
    transmute(
      Factor = factor_label_for_table(condition),
      EQ_Profile = eq,
      Volume_Shift_dB = mean,
      SD_Random_Intercepts_dB = NA_real_,
      SD = sd,
      CI_2_5 = q2_5,
      CI_97_5 = q97_5
    ),
  summarise_plot_stats(volume_random_draws) %>%
    transmute(
      Factor = factor_label_for_table(condition),
      EQ_Profile = eq,
      Volume_Shift_dB = NA_real_,
      SD_Random_Intercepts_dB = mean,
      SD = sd,
      CI_2_5 = q2_5,
      CI_97_5 = q97_5
    )
) %>%
  arrange(match(Factor, factor_order), EQ_Profile) %>%
  mutate(across(c(Volume_Shift_dB, SD_Random_Intercepts_dB, SD, CI_2_5, CI_97_5), ~ round(.x, 2)))

write.csv(eq_table, OUT_EQ_TABLE, row.names = FALSE)
write.csv(volume_table, OUT_VOLUME_TABLE, row.names = FALSE)
cat("Tables exported successfully.\n")


# 6. Final Visualization -------------------------------------------------------
if (RUN_VISUALIZATION) {
  cat("Generating multi-panel visualization...\n")
  
  EPISODE_ORDER <- c("Enjoyment", "Distraction", "Relaxation", "Focus-Motivation")
  NOISE_ORDER <- c("DIESEL")
  EPISODE_LABELS <- c("Enjoyment" = "Enjoyment\n(vs. Null)", "Distraction" = "Distraction\n(vs. Null)", "Relaxation" = "Relaxation   \n(vs. Null)   ", "Focus-Motivation" = "   Focus-Motivation\n   (vs. Null)")
  NOISE_LABELS <- c("DIESEL" = "Diesel\n(vs. EV)")
  EQ_COLOR_MAP <- c("BassCut" = "#fc8d62", "HighCut" = "#8da0cb", "VocalReduce" = "#e78ac3", "V-Shape" = "#a6d854", "VShape" = "#a6d854", "Volume" = "#4d4d4d")
  
  eq_breaks <- seq(0, 9, by = 1)
  eq_breaks_or1 <- eq_breaks[eq_breaks >= 1]
  eq_labels <- function(x) ifelse(abs(x - round(x)) < 1e-8, sprintf("%4.1f", x), "")
  
  volume_breaks <- seq(-1.5, 4.5, by = 0.5)
  volume_breaks_nonneg <- volume_breaks[volume_breaks >= 0]
  volume_labels <- function(x) sprintf("%4.1f", x)
  
  p1 <- create_effect_panel(eq_fixed_draws %>% filter(condition %in% EPISODE_ORDER), eq_fixed_summary %>% filter(condition %in% EPISODE_ORDER), EPISODE_ORDER, c("Enjoyment" = "", "Distraction" = "", "Relaxation" = "", "Focus-Motivation" = ""), c(0, 9), eq_breaks, eq_labels, "Odds Ratio of\nSelection (vs. Flat)", "(a) Listening Episodes\n", EQ_COLOR_MAP, reference_line = 1, show_x_axis_labels = FALSE, show_legend = TRUE)
  p2 <- create_effect_panel(eq_fixed_draws %>% filter(condition %in% NOISE_ORDER), eq_fixed_summary %>% filter(condition %in% NOISE_ORDER), NOISE_ORDER, c("DIESEL" = ""), c(0, 9), eq_breaks, eq_labels, NULL, "(b) Car Cabin\nNoise", EQ_COLOR_MAP, reference_line = 1, show_y_axis_title = FALSE, show_x_axis_labels = FALSE)
  p3 <- create_effect_panel(eq_random_draws %>% filter(condition == "Excerpt" & or >= 1), eq_random_summary %>% filter(condition == "Excerpt") %>% mutate(lower = pmax(lower, 1)), c("Excerpt"), c("Excerpt" = ""), c(0, 9), eq_breaks_or1, eq_labels, "Median Odds Ratio", "(c) Music\nExcerpt", EQ_COLOR_MAP, reference_line = 1, is_random_effect = TRUE, show_x_axis_labels = FALSE)
  p4 <- create_effect_panel(eq_random_draws %>% filter(condition == "Individual" & or >= 1), eq_random_summary %>% filter(condition == "Individual") %>% mutate(lower = pmax(lower, 1)), c("Individual"), c("Individual" = ""), c(0, 9), eq_breaks_or1, eq_labels, NULL, "(d) Participant\n", EQ_COLOR_MAP, reference_line = 1, is_random_effect = TRUE, show_y_axis_title = FALSE, show_x_axis_labels = FALSE)
  
  p5 <- create_effect_panel(volume_fixed_draws %>% filter(condition %in% EPISODE_ORDER), volume_fixed_summary %>% filter(condition %in% EPISODE_ORDER), EPISODE_ORDER, EPISODE_LABELS, c(-1.5, 4.5), volume_breaks, volume_labels, "Estimated Volume\nShift (dB)", NULL, EQ_COLOR_MAP, reference_line = 0)
  p6 <- create_effect_panel(volume_fixed_draws %>% filter(condition %in% NOISE_ORDER), volume_fixed_summary %>% filter(condition %in% NOISE_ORDER), NOISE_ORDER, NOISE_LABELS, c(-1.5, 4.5), volume_breaks, volume_labels, NULL, NULL, EQ_COLOR_MAP, reference_line = 0, show_y_axis_title = FALSE)
  p7 <- create_effect_panel(volume_random_draws %>% filter(condition == "Excerpt" & or >= 0), volume_random_summary %>% filter(condition == "Excerpt") %>% mutate(lower = pmax(lower, 0)), c("Excerpt"), c("Excerpt" = "Music Excerpt\n"), c(-1.5, 4.5), volume_breaks_nonneg, volume_labels, "Standard Deviation of\nRandom Intercepts (dB)", NULL, EQ_COLOR_MAP, reference_line = 0, is_random_effect = TRUE)
  p8 <- create_effect_panel(volume_random_draws %>% filter(condition == "Individual" & or >= 0), volume_random_summary %>% filter(condition == "Individual") %>% mutate(lower = pmax(lower, 0)), c("Individual"), c("Individual" = "Participant\n"), c(-1.5, 4.5), volume_breaks_nonneg, volume_labels, NULL, NULL, EQ_COLOR_MAP, reference_line = 0, is_random_effect = TRUE, show_y_axis_title = FALSE)
  
  pdf(OUT_COMBINED, width = 16, height = 9)
  grid::grid.newpage()
  lay <- grid::grid.layout(nrow = 5, ncol = 9, widths = grid::unit.c(grid::unit(0.15, "in"), grid::unit(3.5, "null"), grid::unit(0.04, "in"), grid::unit(0.8, "null"), grid::unit(0.04, "in"), grid::unit(0.8, "null"), grid::unit(0.04, "in"), grid::unit(0.8, "null"), grid::unit(0.15, "in")), heights = grid::unit.c(grid::unit(0.18, "in"), grid::unit(0, "in"), grid::unit(1, "null"), grid::unit(1, "null"), grid::unit(0.18, "in")))
  grid::pushViewport(grid::viewport(layout = lay))
  
  g1 <- lapply(list(p1, p2, p3, p4), ggplotGrob)
  g2 <- lapply(list(p5, p6, p7, p8), ggplotGrob)
  for (idx in seq_along(g1)) {
    mw <- grid::unit.pmax(g1[[idx]]$widths, g2[[idx]]$widths)
    g1[[idx]]$widths <- mw
    g2[[idx]]$widths <- mw
  }
  
  cols <- c(2, 4, 6, 8)
  for (idx in seq_along(cols)) {
    grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = cols[idx]))
    grid::grid.draw(g1[[idx]])
    grid::upViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 4, layout.pos.col = cols[idx]))
    grid::grid.draw(g2[[idx]])
    grid::upViewport()
  }
  grid::upViewport()
  dev.off()
  cat("Combined multi-panel plot generated successfully.\n")
}
