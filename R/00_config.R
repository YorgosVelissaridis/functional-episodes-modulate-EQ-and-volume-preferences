# ==============================================================================
# 00_config.R: Project Configuration and Dependencies
# ==============================================================================

suppressPackageStartupMessages(library(brms))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(posterior))
suppressPackageStartupMessages(library(ggdist))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(loo))
suppressPackageStartupMessages(library(cmdstanr))

# ------------------------------------------------------------------------------
# 1. GLMM / Bayesian Modeling Configuration
# ------------------------------------------------------------------------------
ITER <- 2000
CHAINS <- 4
WARMUP <- floor(ITER / 2)
ADAPT_DELTA <- 0.90
LOO_CORES <- 1

EQ_FORMULA <- selected_eq ~ noise + episode + (1 | participant_id) + (1 | excerpt_id)
VOLUME_FORMULA <- selected_volume_db ~ noise + episode + (1 | participant_id) + (1 | excerpt_id)

# Confound Models
EQ_FORMULA_FULL <- selected_eq ~ selected_volume_db + noise + episode + (1 | participant_id) + (1 | excerpt_id)

# ------------------------------------------------------------------------------
# 2. LOO / LOCO Configuration
# ------------------------------------------------------------------------------
USE_MOMENT_MATCH <- TRUE
RELOO <- FALSE
AUTO_RELOO_IF_HIGH_K <- FALSE
FUTURE_GLOBALS_MAXSIZE <- 2 * 1024^3

# ------------------------------------------------------------------------------
# 3. Paths and Directories
# ------------------------------------------------------------------------------
DATA_FILE <- "./data/data_en_screened.csv"
RESULTS_DIR <- "./results"
CACHE_DIR <- file.path(RESULTS_DIR, "cache")
MODEL_CACHE_DIR <- file.path(CACHE_DIR, "models")
LOO_CACHE_DIR <- file.path(CACHE_DIR, "loo")

# Main Models
EQ_MODEL_FILE <- file.path(MODEL_CACHE_DIR, "fullmodel_eq.rds")
VOLUME_MODEL_FILE <- file.path(MODEL_CACHE_DIR, "fullmodel_volume.rds")
NULL_EQ_MODEL_FILE <- file.path(MODEL_CACHE_DIR, "nullmodel_eq.rds")

# LOCO Caches
DELTA_CACHE_CSV <- file.path(CACHE_DIR, "elpd_cache.csv") 
DELTA_CACHE_RDS <- file.path(CACHE_DIR, "elpd_cache.rds")

# Outputs
OUT_COMBINED <- file.path(RESULTS_DIR, "fullmodel_plot.pdf")
OUT_ELPD_PLOT <- file.path(RESULTS_DIR, "elpd_plot.pdf")
OUT_EQ_TABLE <- file.path(RESULTS_DIR, "fullmodel_eq_table.csv")
OUT_VOLUME_TABLE <- file.path(RESULTS_DIR, "fullmodel_volume_table.csv")
OUT_R2_TABLE <- file.path(RESULTS_DIR, "fullmodel_r2_table.csv")
OUT_ELPD_TABLE <- file.path(RESULTS_DIR, "loco_elpd_table.csv")

# Ensure output structures exist
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MODEL_CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOO_CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 4. Execution Flags
# ------------------------------------------------------------------------------
RUN_MODELING <- TRUE
RUN_VISUALIZATION <- TRUE
RUN_LOCO_PLOTTING <- TRUE

# Script 02: Main Full Models
FORCE_REFIT_FULL_MODELS <- FALSE
FORCE_RECOMPUTE_FULL_LOO <- FALSE

# Script 02: LOCO Models
FORCE_REFIT_LOCO_MODELS <- FALSE
FORCE_RECOMPUTE_LOCO_LOO <- FALSE
FORCE_RECOMPUTE_DELTA <- FALSE

# Script 04: Confound Models with Auxiliary Covariates
FORCE_REFIT_CONFOUND_MODELS <- FALSE
FORCE_RECOMPUTE_CONFOUND_LOO <- FALSE

# ------------------------------------------------------------------------------
# 5. LOCO Reduced Model Definitions
# ------------------------------------------------------------------------------
MODEL_SPECS <- list(
  list(
    tag = "drop_episode",
    eq = selected_eq ~ noise + (1 | participant_id) + (1 | excerpt_id),
    volume = selected_volume_db ~ noise + (1 | participant_id) + (1 | excerpt_id)
  ),
  list(
    tag = "drop_noise",
    eq = selected_eq ~ episode + (1 | participant_id) + (1 | excerpt_id),
    volume = selected_volume_db ~ episode + (1 | participant_id) + (1 | excerpt_id)
  ),
  list(
    tag = "drop_excerpt",
    eq = selected_eq ~ noise + episode + (1 | participant_id),
    volume = selected_volume_db ~ noise + episode + (1 | participant_id)
  ),
  list(
    tag = "drop_participant",
    eq = selected_eq ~ noise + episode + (1 | excerpt_id),
    volume = selected_volume_db ~ noise + episode + (1 | excerpt_id)
  ),
  list(
    tag = "drop_excerpt_add_octband",
    eq = selected_eq ~ noise + episode +
      oct_level_db_31_5 + oct_level_db_63 + oct_level_db_125 + oct_level_db_250 +
      oct_level_db_500 + oct_level_db_1000 + oct_level_db_2000 + oct_level_db_4000 +
      oct_level_db_8000 + oct_level_db_16000 +
      (1 | participant_id),
    volume = selected_volume_db ~ noise + episode +
      oct_level_db_31_5 + oct_level_db_63 + oct_level_db_125 + oct_level_db_250 +
      oct_level_db_500 + oct_level_db_1000 + oct_level_db_2000 + oct_level_db_4000 +
      oct_level_db_8000 + oct_level_db_16000 +
      (1 | participant_id)
  ),
  list(
    tag = "drop_excerpt_add_vocalpresence",
    eq = selected_eq ~ noise + episode + vocal_presence + (1 | participant_id),
    volume = selected_volume_db ~ noise + episode + vocal_presence + (1 | participant_id)
  ),
  list(
    tag = "drop_excerpt_add_genres",
    eq = selected_eq ~ noise + episode + genres + (1 | participant_id),
    volume = selected_volume_db ~ noise + episode + genres + (1 | participant_id)
  ),
  list(
    tag = "drop_participant_add_age",
    eq = selected_eq ~ noise + episode + age + (1 | excerpt_id),
    volume = selected_volume_db ~ noise + episode + age + (1 | excerpt_id)
  ),
  list(
    tag = "drop_participant_add_gender",
    eq = selected_eq ~ noise + episode + gender_merged + (1 | excerpt_id),
    volume = selected_volume_db ~ noise + episode + gender_merged + (1 | excerpt_id)
  ),
  list(
    tag = "drop_participant_add_personality",
    eq = selected_eq ~ noise + episode +
      Openness + Neuroticism + Extraversion + Conscientiousness + Agreeableness +
      (1 | excerpt_id),
    volume = selected_volume_db ~ noise + episode +
      Openness + Neuroticism + Extraversion + Conscientiousness + Agreeableness +
      (1 | excerpt_id)
  )
)

PLOT_GROUPS <- list(
  list(
    panel = "(a) Impact of Factor Removal\non EQ Prediction",
    model = c("drop_episode", "drop_noise", "drop_excerpt", "drop_participant")
  ),
  list(
    panel = "(b) Predictive Contribution of Acoustic Features\non EQ Prediction",
    model = c("drop_excerpt", "drop_excerpt_add_octband", "drop_excerpt_add_vocalpresence", "drop_excerpt_add_genres")
  ),
  list(
    panel = "(c) Predictive Contribution of Participant Traits\non EQ Prediction",
    model = c("drop_participant", "drop_participant_add_age", "drop_participant_add_gender", "drop_participant_add_personality")
  )
)

MODEL_X_LABELS <- c(
  drop_episode = "w/o\nFunctional\nEpisode",
  drop_noise = "w/o\nCabin\nNoise",
  drop_excerpt = "w/o\nMusic\nExcerpt",
  drop_participant = "w/o\nParticipant",
  drop_excerpt_add_octband = "Replaced\nw/ Oct Band\nEnergy",
  drop_excerpt_add_vocalpresence = "Replaced\nw/ Vocal\nPresence",
  drop_excerpt_add_genres = "Replaced w/\nGenre\nTags",
  drop_participant_add_age = "Replaced\nw/\nAge",
  drop_participant_add_gender = "Replaced\nw/\nGender",
  drop_participant_add_personality = "Replaced\nw/\nPersonality"
)

# ------------------------------------------------------------------------------
# 6. Plot Settings
# ------------------------------------------------------------------------------
AXIS_TICK_TEXT_SIZE <- 13
AXIS_TICK_TEXT_COLOR <- "grey30"
PLOT_TITLE_SIZE <- 14
SUBPLOT_COL_SPACER_IN <- 0.22
Y_AXIS_LABEL_EXPR <- expression(
  atop(
    "Change in Predictive Accuracy",
    Delta * "ELPD (Model - Full Model)"
  )
)