# =============================================================================
# Learning Debt simulation
# Fast, restartable deployment-calibrated simulation with expanded baselines and score-unit sensitivity
#
# This script is meant to be run with Rscript, not inside an interactive
# RStudio session. It checkpoints the expensive stages and can resume after an
# interrupted session.
#
# Main design:
#   1. Warm-start deployed and shadow NIG posteriors from the same initial data.
#   2. Use separate update, monitor, and evaluation batches each period.
#   3. Apply period-t actions after period-t evaluation, so retraining affects
#      period t+1 and later.
#   4. Evaluate policies by lambda * retrains plus accumulated positive
#      evaluation-score regret relative to the shadow posterior.
#   5. Age-adjust sqrt(KL) using no-shift calibration paths.
#   6. Fit debt-only and hybrid utility calibrators for expected predictive
#      regret. Retrain when calibrated expected regret exceeds retraining cost.
#   7. Tune every policy on calibration paths and evaluate on held-out paths.
#
# Environment variables:
#   LD_RUN_MODE    smoke | paper | full     Default: paper
#   LD_OUTPUT_DIR  output directory         Default: ~/
#   LD_RESUME      1 or 0                   Default: 1
#   LD_N_CAL       override calibration paths per scenario cell
#   LD_N_TUNE      override tuning paths per scenario cell
#   LD_N_TEST      override test paths per scenario cell
#
# =============================================================================

options(stringsAsFactors = FALSE)

# =============================================================================
# 0. Configuration
# =============================================================================

RUN_MODE <- tolower(Sys.getenv("LD_RUN_MODE", unset = "paper"))
default_output <- file.path(path.expand("~"), paste0("learning_debt_sim_outputs_v13d_", RUN_MODE))
OUTPUT_DIR <- path.expand(Sys.getenv("LD_OUTPUT_DIR", unset = default_output))
RESUME <- Sys.getenv("LD_RESUME", unset = "1") != "0"

T_PERIODS <- 200L
BURN_IN   <- 20L
N_INIT    <- 300L
N_UPDATE  <- 5L
N_MONITOR <- 5L
N_EVAL    <- 40L

KAPPAS      <- c(0.1, 0.25, 0.5, 1.0, 2.0, 4.0)
SHIFT_PROBS <- c(0.02, 0.05, 0.10, 0.20)
SHIFT_TYPES <- c("none", "abrupt_coef", "variance", "gradual")

if (RUN_MODE == "smoke") {
  T_PERIODS <- 50L
  BURN_IN   <- 8L
  N_INIT    <- 80L
  N_EVAL    <- 30L
  KAPPAS <- c(0.5, 1.0, 2.0)
  SHIFT_PROBS <- c(0.05, 0.10)
  N_CAL_PATHS_PER_CELL  <- 2L
  N_TUNE_PATHS_PER_CELL <- 2L
  N_TEST_PATHS_PER_CELL <- 3L
} else if (RUN_MODE == "full") {
  N_CAL_PATHS_PER_CELL  <- 60L
  N_TUNE_PATHS_PER_CELL <- 25L
  N_TEST_PATHS_PER_CELL <- 300L
} else {
  # Default paper mode. This is the practical submission run for managed
  # RStudio environments. It is large enough to be useful but much less likely
  # to be interrupted than a full 300-rep stress run.
  N_CAL_PATHS_PER_CELL  <- 20L
  N_TUNE_PATHS_PER_CELL <- 8L
  N_TEST_PATHS_PER_CELL <- 100L
}

# Optional overrides.
get_env_int <- function(name, default) {
  x <- Sys.getenv(name, unset = "")
  if (nchar(x) == 0) return(default)
  y <- suppressWarnings(as.integer(x))
  if (!is.finite(y) || y <= 0L) default else y
}
N_CAL_PATHS_PER_CELL  <- get_env_int("LD_N_CAL",  N_CAL_PATHS_PER_CELL)
N_TUNE_PATHS_PER_CELL <- get_env_int("LD_N_TUNE", N_TUNE_PATHS_PER_CELL)
N_TEST_PATHS_PER_CELL <- get_env_int("LD_N_TEST", N_TEST_PATHS_PER_CELL)

# Expanded candidate grids. These are deliberately wider than v13c because the
CAL_INTERVALS        <- c(1L, 2L, 3L, 5L, 10L, 20L, 40L, 80L)
CUSUM_THRESHOLDS     <- c(0.01, 0.025, 0.05, 0.10, 0.20, 0.50, 1, 2, 4, 8, 16)
ALARM_QUANTILES      <- c(0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.975)
UTILITY_SCALES       <- c(0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 4.00, 8.00)
DEBT_THRESHOLD_PROBS <- unique(c(seq(0.50, 0.95, by = 0.05), 0.975, 0.99))

CAL_RANDOM_RESET_PROB <- 1 / 25
SUPPRESS_FINAL_ACTION <- TRUE

SEED_CAL_BASE  <- 120000L
SEED_TUNE_BASE <- 220000L
SEED_TEST_BASE <- 320000L
SEED_ILLUS_BASE <- 440000L

PRIOR <- list(mu_0 = 0, kappa_0 = 1, alpha_0 = 2, beta_0 = 1)

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) {
  cat(sprintf(...), "\n")
  flush.console()
}

# =============================================================================
# 1. NIG posterior algebra
# =============================================================================

bayes_update_xy <- function(prior, x, y) {
  if (length(x) == 0L) return(prior)
  ss_xy <- sum(x * y)
  ss_xx <- sum(x * x)
  n <- length(x)
  kappa_n <- prior$kappa_0 + ss_xx
  mu_n <- (prior$kappa_0 * prior$mu_0 + ss_xy) / kappa_n
  alpha_n <- prior$alpha_0 + n / 2
  beta_n <- prior$beta_0 +
    0.5 * (sum(y * y) - kappa_n * mu_n^2 + prior$kappa_0 * prior$mu_0^2)
  list(mu_0 = mu_n, kappa_0 = kappa_n, alpha_0 = alpha_n, beta_0 = beta_n)
}

kl_nig_exact <- function(p, q) {
  ap <- p$alpha_0; bp <- p$beta_0; kp <- p$kappa_0; mp <- p$mu_0
  aq <- q$alpha_0; bq <- q$beta_0; kq <- q$kappa_0; mq <- q$mu_0
  if (ap <= 0 || bp <= 0 || kp <= 0 || aq <= 0 || bq <= 0 || kq <= 0) return(NA_real_)
  kl_ig <- aq * log(bp / bq) + lgamma(aq) - lgamma(ap) +
    (ap - aq) * digamma(ap) + ap * (bq / bp - 1)
  e_inv_sigma2_p <- ap / bp
  kl_beta <- 0.5 * (log(kp / kq) + kq / kp - 1 + kq * (mp - mq)^2 * e_inv_sigma2_p)
  max(kl_ig + kl_beta, 0)
}

pred_log_score_xy <- function(post, x, y) {
  pred_mean <- post$mu_0 * x
  pred_var <- (post$beta_0 / post$alpha_0) * (1 + x * x / post$kappa_0)
  pred_var <- pmax(pred_var, 1e-12)
  df <- 2 * post$alpha_0
  mean(dt((y - pred_mean) / sqrt(pred_var), df = df, log = TRUE) - log(sqrt(pred_var)))
}

# =============================================================================
# 2. Helpers
# =============================================================================

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, names = FALSE, type = 7))
}

se_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

safe_ratio <- function(num, den) {
  ifelse(is.finite(num) & is.finite(den) & den > 0, num / den, NA_real_)
}

finite_scalar <- function(x, default = 0) {
  if (length(x) == 0L) return(default)
  y <- suppressWarnings(as.numeric(x[1]))
  if (!is.finite(y)) default else y
}

finite_vector <- function(x, default = 0) {
  y <- suppressWarnings(as.numeric(x))
  names(y) <- names(x)
  y[!is.finite(y)] <- default
  y
}

safe_gt <- function(a, b) {
  # Allows Inf thresholds, but never returns NA.
  aa <- suppressWarnings(as.numeric(a[1]))
  bb <- suppressWarnings(as.numeric(b[1]))
  if (length(aa) == 0L || length(bb) == 0L || is.na(aa) || is.na(bb)) return(FALSE)
  aa > bb
}

safe_action <- function(x) {
  if (length(x) == 0L) return(0L)
  y <- suppressWarnings(as.numeric(x[1]))
  if (!is.finite(y) || is.na(y)) return(0L)
  as.integer(y != 0)
}

make_seed <- function(base, shift_type, shift_prob, sim_id) {
  st_id <- match(shift_type, SHIFT_TYPES)
  sp_id <- match(shift_prob, SHIFT_PROBS)
  as.integer(base + 100000L * st_id + 10000L * sp_id + sim_id)
}

append_csv <- function(df, file) {
  write.table(
    df,
    file = file,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(file),
    append = file.exists(file),
    qmethod = "double"
  )
}

checkpoint_file <- function(name) file.path(OUTPUT_DIR, name)

# =============================================================================
# 3. Data generation
# =============================================================================

generate_path <- function(T, shift_prob, shift_type) {
  G_true <- integer(T)
  beta_vec <- numeric(T)
  sigma2_vec <- rep(1, T)
  cur_beta <- 0
  cur_sigma2 <- 1

  for (t in seq_len(T)) {
    prev_G <- if (t == 1L) 0L else G_true[t - 1L]

    if (shift_type == "none") {
      G_true[t] <- 0L
    } else if (shift_type == "abrupt_coef") {
      if (t > BURN_IN && prev_G == 0L && runif(1) < shift_prob) {
        G_true[t] <- 1L
        cur_beta <- cur_beta + rnorm(1, 0, 2)
      } else {
        G_true[t] <- prev_G
      }
    } else if (shift_type == "variance") {
      if (t > BURN_IN && prev_G == 0L && runif(1) < shift_prob) {
        G_true[t] <- 1L
        cur_sigma2 <- cur_sigma2 * runif(1, 3, 6)
      } else {
        G_true[t] <- prev_G
      }
    } else if (shift_type == "gradual") {
      if (t > BURN_IN && prev_G == 0L && runif(1) < shift_prob) {
        G_true[t] <- 1L
      } else {
        G_true[t] <- prev_G
      }
      if (G_true[t] == 1L) cur_beta <- cur_beta + rnorm(1, 0, 0.15)
    } else {
      stop("Unknown shift_type: ", shift_type)
    }

    beta_vec[t] <- cur_beta
    sigma2_vec[t] <- cur_sigma2
  }

  make_matrix <- function(n) matrix(rnorm(T * n), nrow = T, ncol = n)
  xu <- make_matrix(N_UPDATE)
  xm <- make_matrix(N_MONITOR)
  xe <- make_matrix(N_EVAL)
  yu <- beta_vec * xu + sqrt(sigma2_vec) * make_matrix(N_UPDATE)
  ym <- beta_vec * xm + sqrt(sigma2_vec) * make_matrix(N_MONITOR)
  ye <- beta_vec * xe + sqrt(sigma2_vec) * make_matrix(N_EVAL)

  x0 <- rnorm(N_INIT)
  y0 <- rnorm(N_INIT)

  list(
    init_x = x0, init_y = y0,
    G_true = G_true, beta = beta_vec, sigma2 = sigma2_vec,
    xu = xu, yu = yu, xm = xm, ym = ym, xe = xe, ye = ye
  )
}

# =============================================================================
# 4. Feature computation and calibration models
# =============================================================================

compute_features <- function(shadow_post, dep_post, path, t, spell_age) {
  d_val <- kl_nig_exact(shadow_post, dep_post)
  if (!is.finite(d_val)) d_val <- 0

  ls_shadow_monitor <- pred_log_score_xy(shadow_post, path$xm[t, ], path$ym[t, ])
  ls_dep_monitor <- pred_log_score_xy(dep_post, path$xm[t, ], path$ym[t, ])
  ls_shadow_eval <- pred_log_score_xy(shadow_post, path$xe[t, ], path$ye[t, ])
  ls_dep_eval <- pred_log_score_xy(dep_post, path$xe[t, ], path$ye[t, ])

  score_gap <- ls_shadow_monitor - ls_dep_monitor
  if (!is.finite(score_gap)) score_gap <- 0
  eval_gap <- ls_shadow_eval - ls_dep_eval
  if (!is.finite(eval_gap)) eval_gap <- 0

  dep_scale <- sqrt(max(dep_post$beta_0 / dep_post$alpha_0, 1e-12))
  resid_exceed <- mean(abs(path$ym[t, ] - dep_post$mu_0 * path$xm[t, ]) > 2 * dep_scale)

  c(
    spell_age = spell_age,
    sqrt_debt = sqrt(max(d_val, 0)),
    debt_exact = d_val,
    score_gap = score_gap,
    pos_score_gap = max(score_gap, 0),
    param_gap = abs(shadow_post$mu_0 - dep_post$mu_0),
    resid_exceed = resid_exceed,
    eval_score_gap = eval_gap,
    regret_eval = max(eval_gap, 0)
  )
}

collect_calibration_rows <- function(n_paths_per_cell) {
  cache <- checkpoint_file("calibration_rows_raw_v13d.rds")
  if (RESUME && file.exists(cache)) {
    log_msg("Loading cached calibration rows: %s", cache)
    return(readRDS(cache))
  }

  grid <- expand.grid(
    shift_type = SHIFT_TYPES,
    shift_prob = SHIFT_PROBS,
    sim_id = seq_len(n_paths_per_cell),
    stringsAsFactors = FALSE
  )

  log_msg("Collecting calibration rows: %d paths", nrow(grid))
  rows <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    if (i %% 25L == 0L) log_msg("  calibration path %d / %d", i, nrow(grid))
    r <- grid[i, ]
    set.seed(make_seed(SEED_CAL_BASE, r$shift_type, r$shift_prob, r$sim_id))
    path <- generate_path(T_PERIODS, r$shift_prob, r$shift_type)

    init_post <- bayes_update_xy(PRIOR, path$init_x, path$init_y)
    dep_post <- init_post
    shadow_post <- init_post
    last_retrain <- 0L

    mat <- matrix(NA_real_, nrow = T_PERIODS, ncol = 9L)
    colnames(mat) <- c("spell_age", "sqrt_debt", "debt_exact", "score_gap",
                       "pos_score_gap", "param_gap", "resid_exceed",
                       "eval_score_gap", "regret_eval")

    for (t in seq_len(T_PERIODS)) {
      shadow_post <- bayes_update_xy(shadow_post, path$xu[t, ], path$yu[t, ])
      spell_age <- t - last_retrain
      mat[t, ] <- compute_features(shadow_post, dep_post, path, t, spell_age)
      if (t < T_PERIODS && runif(1) < CAL_RANDOM_RESET_PROB) {
        dep_post <- shadow_post
        last_retrain <- t
      }
    }

    rows[[i]] <- data.frame(
      shift_type = r$shift_type,
      shift_prob = r$shift_prob,
      sim_id = r$sim_id,
      t = seq_len(T_PERIODS),
      G_true = path$G_true,
      beta_t = path$beta,
      sigma2_t = path$sigma2,
      mat,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  saveRDS(out, cache)
  out
}

fit_age_adjuster <- function(cal_df) {
  stable_df <- cal_df[cal_df$shift_type == "none", ]
  if (nrow(stable_df) < 50L) stable_df <- cal_df[cal_df$G_true == 0L, ]
  if (nrow(stable_df) < 50L) stable_df <- cal_df
  stable_df$log_age <- log1p(stable_df$spell_age)
  stable_df$sqrt_age <- sqrt(stable_df$spell_age)
  m <- lm(sqrt_debt ~ log_age + sqrt_age + spell_age, data = stable_df)
  s <- sd(residuals(m), na.rm = TRUE)
  if (!is.finite(s) || s < 0.01) s <- 0.01
  list(coef = coef(m), resid_sd = s, model = m)
}

manual_age_pred <- function(age_adjuster, spell_age) {
  co <- age_adjuster$coef
  getc <- function(nm) {
    if (nm %in% names(co) && is.finite(co[[nm]])) co[[nm]] else 0
  }
  out <- getc("(Intercept)") + getc("log_age") * log1p(spell_age) +
    getc("sqrt_age") * sqrt(spell_age) + getc("spell_age") * spell_age
  as.numeric(out)
}

add_calibrated_features <- function(df, age_adjuster) {
  df$log_age <- log1p(df$spell_age)
  df$sqrt_age <- sqrt(df$spell_age)
  df$debt_age_mean <- manual_age_pred(age_adjuster, df$spell_age)
  df$debt_z <- (df$sqrt_debt - df$debt_age_mean) / age_adjuster$resid_sd
  df$debt_z_pos <- pmax(df$debt_z, 0)
  df$debt_z_pos2 <- df$debt_z_pos^2
  df$log_age2 <- df$log_age^2
  df$log_regret <- log1p(pmax(df$regret_eval, 0))
  df$pos_score_gap <- pmax(df$score_gap, 0)
  df
}

fit_utility_models <- function(cal_df_adj) {
  model_df <- cal_df_adj[is.finite(cal_df_adj$log_regret) & is.finite(cal_df_adj$debt_z_pos), ]
  if (nrow(model_df) < 100L) stop("Too few calibration rows to fit utility models.")

  debt_m <- lm(log_regret ~ debt_z_pos + debt_z_pos2 + log_age + log_age2, data = model_df)
  hybrid_m <- lm(log_regret ~ debt_z_pos + debt_z_pos2 + log_age + log_age2 +
                   pos_score_gap + param_gap + resid_exceed + sqrt_debt,
                 data = model_df)
  cap <- safe_quantile(model_df$regret_eval, 0.995)
  if (!is.finite(cap) || cap <= 0) cap <- max(model_df$regret_eval, na.rm = TRUE)
  if (!is.finite(cap) || cap <= 0) cap <- 1
  list(
    debt_coef = coef(debt_m),
    hybrid_coef = coef(hybrid_m),
    debt_model = debt_m,
    hybrid_model = hybrid_m,
    regret_cap = cap
  )
}

manual_lm_predict_one <- function(co, x) {
  s <- 0
  for (nm in names(co)) {
    b <- co[[nm]]
    if (!is.finite(b)) next
    if (nm == "(Intercept)") {
      s <- s + b
    } else if (nm %in% names(x)) {
      xv <- suppressWarnings(as.numeric(x[[nm]][1]))
      if (!is.finite(xv)) next
      s <- s + b * xv
    }
  }
  if (!is.finite(s)) 0 else s
}

predict_expected_regret_one <- function(model_bundle, feat, policy) {
  co <- if (policy == "debt_utility") model_bundle$debt_coef else model_bundle$hybrid_coef
  lp <- manual_lm_predict_one(co, feat)
  if (!is.finite(lp)) lp <- 0
  cap <- finite_scalar(model_bundle$regret_cap, 1)
  if (cap <= 0) cap <- 1
  raw <- expm1(lp)
  if (!is.finite(raw)) raw <- cap
  min(max(raw, 0), cap)
}

make_feature_adj_one <- function(feat_vec, age_adjuster) {
  # Named numeric vector in, named numeric vector out. All nonfinite quantities
  # are made harmless here so a rare numerical irregularity cannot produce an
  # NA policy action during tuning.
  feat_vec <- finite_vector(feat_vec, default = 0)
  spell_age <- max(finite_scalar(feat_vec["spell_age"], 0), 0)
  log_age <- log1p(spell_age)
  debt_age_mean <- finite_scalar(manual_age_pred(age_adjuster, spell_age), 0)
  resid_sd <- finite_scalar(age_adjuster$resid_sd, 0.01)
  if (resid_sd <= 0) resid_sd <- 0.01
  debt_z <- (finite_scalar(feat_vec["sqrt_debt"], 0) - debt_age_mean) / resid_sd
  if (!is.finite(debt_z)) debt_z <- 0
  debt_z_pos <- max(debt_z, 0)
  out <- c(
    feat_vec,
    log_age = log_age,
    sqrt_age = sqrt(spell_age),
    debt_age_mean = debt_age_mean,
    debt_z = debt_z,
    debt_z_pos = debt_z_pos,
    debt_z_pos2 = debt_z_pos^2,
    log_age2 = log_age^2,
    log_regret = log1p(max(finite_scalar(feat_vec["regret_eval"], 0), 0))
  )
  finite_vector(out, default = 0)
}

# =============================================================================
# 5. Policy simulation
# =============================================================================

simulate_policy <- function(path, policy, lambda, age_adjuster, utility_models,
                            param_value = NA_real_, return_series = FALSE) {
  policy <- as.character(policy[1])
  lambda <- finite_scalar(lambda, 1)
  param_value <- finite_scalar(param_value, NA_real_)

  init_post <- bayes_update_xy(PRIOR, path$init_x, path$init_y)
  dep_post <- init_post
  shadow_post <- init_post
  last_retrain <- 0L

  total_regret <- 0
  n_retrain <- 0L
  cusum_val <- 0
  alarm_history <- numeric(0)

  if (return_series) series <- vector("list", T_PERIODS)

  for (t in seq_len(T_PERIODS)) {
    shadow_post <- bayes_update_xy(shadow_post, path$xu[t, ], path$yu[t, ])
    spell_age <- t - last_retrain
    feat <- compute_features(shadow_post, dep_post, path, t, spell_age)
    feat_adj <- make_feature_adj_one(feat, age_adjuster)

    total_regret <- total_regret + finite_scalar(feat_adj["regret_eval"], 0)
    action <- 0L

    if (policy == "calendar") {
      interval <- as.integer(round(finite_scalar(param_value, 10)))
      if (!is.finite(interval) || interval <= 0L) interval <- 10L
      action <- safe_action(t %% interval == 0L)
    } else if (policy == "cusum") {
      threshold <- finite_scalar(param_value, 4)
      if (!is.finite(threshold) || threshold <= 0) threshold <- 4
      cusum_val <- max(0, finite_scalar(cusum_val, 0) + finite_scalar(feat_adj["pos_score_gap"], 0))
      action <- safe_action(safe_gt(cusum_val, threshold))
    } else if (policy == "alarm") {
      q <- finite_scalar(param_value, 0.90)
      if (!is.finite(q) || q <= 0 || q >= 1) q <- 0.90
      alarm_thresh <- if (length(alarm_history) >= 10L) safe_quantile(alarm_history, q) else Inf
      action <- safe_action(safe_gt(feat_adj["debt_z"], alarm_thresh))
      alarm_history <- c(alarm_history, finite_scalar(feat_adj["debt_z"], 0))
    } else if (policy == "debt_threshold") {
      thresh <- finite_scalar(param_value, 1)
      action <- safe_action(safe_gt(feat_adj["debt_z"], thresh))
    } else if (policy == "always_retrain") {
      action <- 1L
    } else if (policy == "never_retrain") {
      action <- 0L
    } else if (policy == "debt_utility" || policy == "hybrid_utility") {
      scale <- finite_scalar(param_value, 1)
      if (!is.finite(scale) || scale <= 0) scale <- 1
      rhat <- finite_scalar(predict_expected_regret_one(utility_models, feat_adj, policy), 0)
      action <- safe_action(safe_gt(scale * rhat, lambda))
    } else {
      stop("Unknown policy: ", policy)
    }

    if (SUPPRESS_FINAL_ACTION && t == T_PERIODS) action <- 0L
    action <- safe_action(action)

    if (isTRUE(action == 1L)) {
      n_retrain <- n_retrain + 1L
      dep_post <- shadow_post
      last_retrain <- t
      cusum_val <- 0
      alarm_history <- numeric(0)
    }

    if (return_series) {
      rhat <- if (policy == "debt_utility" || policy == "hybrid_utility") {
        predict_expected_regret_one(utility_models, feat_adj, policy)
      } else {
        NA_real_
      }
      series[[t]] <- data.frame(
        t = t,
        G_true = path$G_true[t],
        beta_t = path$beta[t],
        sigma2_t = path$sigma2[t],
        action = action,
        lambda = lambda,
        rhat = rhat,
        t(as.matrix(feat_adj)),
        stringsAsFactors = FALSE
      )
    }
  }

  obj <- finite_scalar(lambda, 0) * n_retrain + finite_scalar(total_regret, 0)
  out <- list(objective = obj, total_regret = finite_scalar(total_regret, 0), n_retrain = n_retrain)
  if (return_series) out$series <- do.call(rbind, series)
  out
}

# =============================================================================
# 6. Tuning

# =============================================================================
# 6. Tuning with expanded baselines and score-unit sensitivity
# =============================================================================

build_score_units <- function(cal_df_adj) {
  positive_regrets <- cal_df_adj$regret_eval[cal_df_adj$regret_eval > 0 & is.finite(cal_df_adj$regret_eval)]
  if (length(positive_regrets) == 0L) positive_regrets <- 1

  med <- stats::median(positive_regrets, na.rm = TRUE)
  q75 <- safe_quantile(positive_regrets, 0.75)
  mn  <- mean(positive_regrets, na.rm = TRUE)

  if (!is.finite(med) || med <= 0) med <- 1
  if (!is.finite(q75) || q75 <= 0) q75 <- med
  if (!is.finite(mn)  || mn  <= 0) mn  <- med

  data.frame(
    score_unit_mode = c("median", "q75", "mean"),
    score_unit = c(med, q75, mn),
    n_positive_regrets = length(positive_regrets),
    median_positive_regret = med,
    mean_positive_regret = mn,
    q25_positive_regret = safe_quantile(positive_regrets, 0.25),
    q75_positive_regret = q75,
    run_mode = RUN_MODE,
    stringsAsFactors = FALSE
  )
}

build_candidate_grid <- function(cal_df_adj) {
  debt_thresholds <- stats::quantile(cal_df_adj$debt_z, probs = DEBT_THRESHOLD_PROBS,
                                     names = FALSE, na.rm = TRUE, type = 7)
  debt_thresholds <- sort(unique(debt_thresholds[is.finite(debt_thresholds)]))
  if (length(debt_thresholds) == 0L) debt_thresholds <- c(0, 1, 2)

  out <- rbind(
    data.frame(policy = "never_retrain", parameter = "none", value = 0),
    data.frame(policy = "always_retrain", parameter = "none", value = 0),
    data.frame(policy = "calendar", parameter = "interval", value = CAL_INTERVALS),
    data.frame(policy = "cusum", parameter = "threshold", value = CUSUM_THRESHOLDS),
    data.frame(policy = "alarm", parameter = "quantile", value = ALARM_QUANTILES),
    data.frame(policy = "debt_threshold", parameter = "threshold", value = debt_thresholds),
    data.frame(policy = "debt_utility", parameter = "scale", value = UTILITY_SCALES),
    data.frame(policy = "hybrid_utility", parameter = "scale", value = UTILITY_SCALES)
  )
  out$candidate_id <- seq_len(nrow(out))
  out
}

lambda_for <- function(score_units_df, mode_index, kappa) {
  finite_scalar(kappa, 1) * finite_scalar(score_units_df$score_unit[mode_index], 1)
}

policy_depends_on_lambda <- function(policy) {
  policy %in% c("debt_utility", "hybrid_utility")
}

reobjective <- function(res, lambda) {
  finite_scalar(lambda, 0) * finite_scalar(res$n_retrain, 0) + finite_scalar(res$total_regret, Inf)
}

tune_policies <- function(candidate_grid, score_units_df, age_adjuster, utility_models,
                          n_paths_per_cell) {
  cache <- checkpoint_file("tuned_policies_v13d.rds")
  state_file <- checkpoint_file("tuning_state_v13d.rds")
  if (RESUME && file.exists(cache)) {
    log_msg("Loading cached tuned policies: %s", cache)
    return(readRDS(cache))
  }

  tune_grid <- expand.grid(
    shift_type = SHIFT_TYPES,
    shift_prob = SHIFT_PROBS,
    sim_id = seq_len(n_paths_per_cell),
    stringsAsFactors = FALSE
  )

  nM <- nrow(score_units_df)
  nK <- length(KAPPAS)
  nC <- nrow(candidate_grid)

  if (RESUME && file.exists(state_file)) {
    st <- readRDS(state_file)
    start_i <- st$next_i
    obj_sum <- st$obj_sum
    regret_sum <- st$regret_sum
    retrain_sum <- st$retrain_sum
    count_arr <- st$count_arr
    log_msg("Resuming tuning at path %d / %d", start_i, nrow(tune_grid))
  } else {
    start_i <- 1L
    obj_sum <- array(0, dim = c(nM, nK, nC))
    regret_sum <- array(0, dim = c(nM, nK, nC))
    retrain_sum <- array(0, dim = c(nM, nK, nC))
    count_arr <- array(0, dim = c(nM, nK, nC))
  }

  log_msg("Tuning policies: %d paths, %d candidates, %d kappas, %d score-unit modes",
          nrow(tune_grid), nC, nK, nM)

  if (start_i <= nrow(tune_grid)) {
    for (i in start_i:nrow(tune_grid)) {
      if (i %% 10L == 0L) log_msg("  tuning path %d / %d", i, nrow(tune_grid))
      r <- tune_grid[i, ]
      set.seed(make_seed(SEED_TUNE_BASE, r$shift_type, r$shift_prob, r$sim_id))
      path <- generate_path(T_PERIODS, r$shift_prob, r$shift_type)

      for (cj in seq_len(nC)) {
        cand <- candidate_grid[cj, ]
        pol <- as.character(cand$policy)
        val <- cand$value

        if (policy_depends_on_lambda(pol)) {
          for (mi in seq_len(nM)) {
            for (ki in seq_along(KAPPAS)) {
              lambda <- lambda_for(score_units_df, mi, KAPPAS[ki])
              res <- tryCatch(
                simulate_policy(path, pol, lambda, age_adjuster, utility_models, val, return_series = FALSE),
                error = function(e) {
                  warning(sprintf("Tuning candidate failed: policy=%s, value=%s, mode=%s, kappa=%s, shift=%s, p=%s, sim=%s: %s",
                                  pol, val, score_units_df$score_unit_mode[mi], KAPPAS[ki],
                                  r$shift_type, r$shift_prob, r$sim_id, e$message))
                  list(objective = Inf, total_regret = Inf, n_retrain = NA_integer_)
                }
              )
              obj_sum[mi, ki, cj] <- obj_sum[mi, ki, cj] + finite_scalar(res$objective, Inf)
              regret_sum[mi, ki, cj] <- regret_sum[mi, ki, cj] + finite_scalar(res$total_regret, Inf)
              retrain_sum[mi, ki, cj] <- retrain_sum[mi, ki, cj] + finite_scalar(res$n_retrain, 0)
              count_arr[mi, ki, cj] <- count_arr[mi, ki, cj] + 1L
            }
          }
        } else {
          # For calendar, CUSUM, alarm, debt-threshold, always, and never, the
          # action path does not depend on lambda. Simulate once and re-objective
          # for every cost scale.
          res0 <- tryCatch(
            simulate_policy(path, pol, 0, age_adjuster, utility_models, val, return_series = FALSE),
            error = function(e) {
              warning(sprintf("Tuning candidate failed: policy=%s, value=%s, shift=%s, p=%s, sim=%s: %s",
                              pol, val, r$shift_type, r$shift_prob, r$sim_id, e$message))
              list(objective = Inf, total_regret = Inf, n_retrain = NA_integer_)
            }
          )
          reg_val <- finite_scalar(res0$total_regret, Inf)
          ret_val <- finite_scalar(res0$n_retrain, 0)
          for (mi in seq_len(nM)) {
            for (ki in seq_along(KAPPAS)) {
              lambda <- lambda_for(score_units_df, mi, KAPPAS[ki])
              obj_sum[mi, ki, cj] <- obj_sum[mi, ki, cj] + lambda * ret_val + reg_val
              regret_sum[mi, ki, cj] <- regret_sum[mi, ki, cj] + reg_val
              retrain_sum[mi, ki, cj] <- retrain_sum[mi, ki, cj] + ret_val
              count_arr[mi, ki, cj] <- count_arr[mi, ki, cj] + 1L
            }
          }
        }
      }

      if (i %% 10L == 0L || i == nrow(tune_grid)) {
        saveRDS(list(next_i = i + 1L, obj_sum = obj_sum, regret_sum = regret_sum,
                     retrain_sum = retrain_sum, count_arr = count_arr), state_file)
      }
    }
  }

  rows <- vector("list", nM * nK * nC)
  idx <- 0L
  for (mi in seq_len(nM)) {
    for (ki in seq_along(KAPPAS)) {
      for (cj in seq_len(nC)) {
        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          score_unit_mode = score_units_df$score_unit_mode[mi],
          score_unit = score_units_df$score_unit[mi],
          kappa = KAPPAS[ki],
          lambda = lambda_for(score_units_df, mi, KAPPAS[ki]),
          candidate_id = candidate_grid$candidate_id[cj],
          policy = candidate_grid$policy[cj],
          parameter = candidate_grid$parameter[cj],
          value = candidate_grid$value[cj],
          mean_objective = obj_sum[mi, ki, cj] / count_arr[mi, ki, cj],
          mean_regret = regret_sum[mi, ki, cj] / count_arr[mi, ki, cj],
          mean_retrains = retrain_sum[mi, ki, cj] / count_arr[mi, ki, cj],
          n_paths = count_arr[mi, ki, cj],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  tuning_summary <- do.call(rbind, rows)

  tuned_list <- list()
  idx <- 0L
  for (mode in score_units_df$score_unit_mode) {
    for (k in KAPPAS) {
      for (pol in unique(candidate_grid$policy)) {
        tmp <- tuning_summary[tuning_summary$score_unit_mode == mode &
                                tuning_summary$kappa == k &
                                tuning_summary$policy == pol, ]
        best <- tmp[which.min(tmp$mean_objective), ]
        idx <- idx + 1L
        tuned_list[[idx]] <- best
      }
    }
  }
  tuned_params <- do.call(rbind, tuned_list)
  tuned_params <- tuned_params[, c("score_unit_mode", "score_unit", "kappa", "lambda", "policy",
                                   "parameter", "value", "mean_objective", "mean_regret",
                                   "mean_retrains", "n_paths")]

  out <- list(tuning_summary = tuning_summary, tuned_params = tuned_params)
  saveRDS(out, cache)
  out
}

get_tuned_value <- function(tuned_params, score_unit_mode, kappa, policy, default_value) {
  row <- tuned_params[tuned_params$score_unit_mode == score_unit_mode &
                        tuned_params$kappa == kappa &
                        tuned_params$policy == policy, ]
  if (nrow(row) == 0L) return(default_value)
  val <- row$value[1]
  if (!is.finite(val)) default_value else val
}

make_boundary_summary <- function(tuned_params, candidate_grid) {
  rows <- vector("list", nrow(tuned_params))
  for (i in seq_len(nrow(tuned_params))) {
    row <- tuned_params[i, ]
    cand <- candidate_grid[candidate_grid$policy == row$policy, ]
    vals <- cand$value[is.finite(cand$value)]
    if (length(vals) <= 1L) {
      min_val <- max_val <- row$value
      lower <- upper <- FALSE
    } else {
      min_val <- min(vals)
      max_val <- max(vals)
      lower <- isTRUE(abs(row$value - min_val) <= 1e-12)
      upper <- isTRUE(abs(row$value - max_val) <= 1e-12)
    }
    rows[[i]] <- data.frame(
      score_unit_mode = row$score_unit_mode,
      kappa = row$kappa,
      policy = row$policy,
      parameter = row$parameter,
      tuned_value = row$value,
      grid_min = min_val,
      grid_max = max_val,
      n_grid_values = length(unique(vals)),
      at_lower_boundary = lower,
      at_upper_boundary = upper,
      at_any_boundary = lower || upper,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# =============================================================================
# 7. Held-out evaluation, checkpointed by path
# =============================================================================

run_test_grid <- function(tuned_params, score_units_df, age_adjuster, utility_models,
                          n_paths_per_cell) {
  all_runs_file <- checkpoint_file("all_runs_v13d.csv")
  completed_file <- checkpoint_file("completed_test_paths_v13d.rds")

  test_grid <- expand.grid(
    shift_type = SHIFT_TYPES,
    shift_prob = SHIFT_PROBS,
    sim_id = seq_len(n_paths_per_cell),
    stringsAsFactors = FALSE
  )
  test_grid$path_key <- paste(test_grid$shift_type, test_grid$shift_prob, test_grid$sim_id, sep = "|")

  if (RESUME && file.exists(completed_file)) {
    completed <- readRDS(completed_file)
    log_msg("Resuming held-out evaluation with %d completed paths", length(completed))
  } else if (RESUME && file.exists(all_runs_file)) {
    old_runs <- read.csv(all_runs_file, stringsAsFactors = FALSE)
    completed <- unique(old_runs$path_key)
    saveRDS(completed, completed_file)
    log_msg("Recovered %d completed paths from existing all_runs file", length(completed))
  } else {
    completed <- character(0)
    if (file.exists(all_runs_file) && !RESUME) file.remove(all_runs_file)
  }

  policies_to_run <- c("never_retrain", "always_retrain", "calendar10", "calendar_tuned",
                       "cusum_tuned", "alarm_tuned", "debt_threshold_tuned",
                       "debt_utility_tuned", "hybrid_utility_tuned")

  log_msg("Held-out evaluation: %d paths, %d score-unit modes, %d kappas, %d policies",
          nrow(test_grid), nrow(score_units_df), length(KAPPAS), length(policies_to_run))

  for (i in seq_len(nrow(test_grid))) {
    r <- test_grid[i, ]
    if (r$path_key %in% completed) next
    if (i %% 20L == 0L) log_msg("  test path %d / %d", i, nrow(test_grid))

    set.seed(make_seed(SEED_TEST_BASE, r$shift_type, r$shift_prob, r$sim_id))
    path <- generate_path(T_PERIODS, r$shift_prob, r$shift_type)

    cache_results <- list()
    simulate_cached <- function(policy, lambda, param, lambda_sensitive) {
      if (lambda_sensitive) {
        key <- paste(policy, format(param, digits = 12), format(lambda, digits = 12), sep = "|")
        run_lambda <- lambda
      } else {
        key <- paste(policy, format(param, digits = 12), "lambda_free", sep = "|")
        run_lambda <- 0
      }
      if (!key %in% names(cache_results)) {
        cache_results[[key]] <<- tryCatch(
          simulate_policy(path, policy, run_lambda, age_adjuster, utility_models, param, FALSE),
          error = function(e) {
            warning(sprintf("Policy failed: policy=%s, lambda=%s, param=%s, shift=%s, p=%s, sim=%s: %s",
                            policy, lambda, param, r$shift_type, r$shift_prob, r$sim_id, e$message))
            list(objective = NA_real_, total_regret = NA_real_, n_retrain = NA_integer_)
          }
        )
      }
      cache_results[[key]]
    }

    path_rows <- list()
    idx <- 0L

    for (mi in seq_len(nrow(score_units_df))) {
      mode <- score_units_df$score_unit_mode[mi]
      score_unit <- score_units_df$score_unit[mi]
      for (kappa in KAPPAS) {
        lambda <- kappa * score_unit
        specs <- list(
          never_retrain = list(policy = "never_retrain", param = 0),
          always_retrain = list(policy = "always_retrain", param = 0),
          calendar10 = list(policy = "calendar", param = 10),
          calendar_tuned = list(policy = "calendar", param = get_tuned_value(tuned_params, mode, kappa, "calendar", 10)),
          cusum_tuned = list(policy = "cusum", param = get_tuned_value(tuned_params, mode, kappa, "cusum", 4)),
          alarm_tuned = list(policy = "alarm", param = get_tuned_value(tuned_params, mode, kappa, "alarm", 0.90)),
          debt_threshold_tuned = list(policy = "debt_threshold", param = get_tuned_value(tuned_params, mode, kappa, "debt_threshold", 1)),
          debt_utility_tuned = list(policy = "debt_utility", param = get_tuned_value(tuned_params, mode, kappa, "debt_utility", 1)),
          hybrid_utility_tuned = list(policy = "hybrid_utility", param = get_tuned_value(tuned_params, mode, kappa, "hybrid_utility", 1))
        )

        for (nm in names(specs)) {
          spec <- specs[[nm]]
          lambda_sensitive <- policy_depends_on_lambda(spec$policy)
          res <- simulate_cached(spec$policy, lambda, spec$param, lambda_sensitive)
          obj <- if (lambda_sensitive) finite_scalar(res$objective, NA_real_) else {
            if (is.finite(finite_scalar(res$total_regret, NA_real_))) reobjective(res, lambda) else NA_real_
          }
          idx <- idx + 1L
          path_rows[[idx]] <- data.frame(
            score_unit_mode = mode,
            score_unit = score_unit,
            kappa = kappa,
            lambda = lambda,
            shift_type = r$shift_type,
            shift_prob = r$shift_prob,
            sim_id = r$sim_id,
            path_key = r$path_key,
            policy = nm,
            base_policy = spec$policy,
            tuned_value = spec$param,
            objective = obj,
            total_regret = res$total_regret,
            n_retrain = res$n_retrain,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    append_csv(do.call(rbind, path_rows), all_runs_file)
    completed <- c(completed, r$path_key)
    if (length(completed) %% 20L == 0L || i == nrow(test_grid)) saveRDS(completed, completed_file)
  }

  read.csv(all_runs_file, stringsAsFactors = FALSE)
}

# =============================================================================
# 8. Summaries
# =============================================================================

summarise_by_group <- function(df, group_cols) {
  key <- do.call(paste, c(df[group_cols], sep = "|"))
  groups <- split(df, key)
  rows <- lapply(groups, function(g) {
    base <- g[1, group_cols, drop = FALSE]
    data.frame(
      base,
      n_rep = sum(is.finite(g$objective)),
      mean_objective = mean(g$objective, na.rm = TRUE),
      se_objective = se_mean(g$objective),
      q025_objective = safe_quantile(g$objective, 0.025),
      median_objective = safe_quantile(g$objective, 0.5),
      q975_objective = safe_quantile(g$objective, 0.975),
      mean_regret = mean(g$total_regret, na.rm = TRUE),
      mean_retrains = mean(g$n_retrain, na.rm = TRUE),
      se_retrains = se_mean(g$n_retrain),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_relative_summary <- function(policy_summary) {
  id_cols <- c("score_unit_mode", "score_unit", "kappa", "lambda", "shift_prob", "shift_type")
  objective_wide <- reshape(policy_summary[, c(id_cols, "policy", "mean_objective")],
                            idvar = id_cols, timevar = "policy", direction = "wide")
  names(objective_wide) <- gsub("mean_objective\\.", "obj_", names(objective_wide))
  regret_wide <- reshape(policy_summary[, c(id_cols, "policy", "mean_regret")],
                         idvar = id_cols, timevar = "policy", direction = "wide")
  names(regret_wide) <- gsub("mean_regret\\.", "regret_", names(regret_wide))
  retrain_wide <- reshape(policy_summary[, c(id_cols, "policy", "mean_retrains")],
                          idvar = id_cols, timevar = "policy", direction = "wide")
  names(retrain_wide) <- gsub("mean_retrains\\.", "retrains_", names(retrain_wide))

  out <- merge(objective_wide, regret_wide, by = id_cols, all = TRUE)
  out <- merge(out, retrain_wide, by = id_cols, all = TRUE)

  add_ratio <- function(target, bm, out_df) {
    target_col <- paste0("obj_", target)
    bm_col <- paste0("obj_", bm)
    rel_col <- paste0("rel_", gsub("_tuned$", "", target), "_vs_", gsub("_tuned$", "", bm))
    if (target_col %in% names(out_df) && bm_col %in% names(out_df)) {
      out_df[[rel_col]] <- safe_ratio(out_df[[target_col]], out_df[[bm_col]])
    }
    out_df
  }

  targets <- c("debt_threshold_tuned", "debt_utility_tuned", "hybrid_utility_tuned")
  benchmarks <- c("calendar10", "calendar_tuned", "cusum_tuned", "alarm_tuned",
                  "debt_threshold_tuned", "always_retrain", "never_retrain")
  for (target in targets) {
    for (bm in benchmarks) {
      if (target != bm) out <- add_ratio(target, bm, out)
    }
  }
  out
}

make_cell_win_summary <- function(relative_summary) {
  targets <- c("debt_threshold_tuned", "debt_utility_tuned", "hybrid_utility_tuned")
  benchmarks <- c("calendar10", "calendar_tuned", "cusum_tuned", "alarm_tuned",
                  "debt_threshold_tuned", "always_retrain", "never_retrain")
  df <- relative_summary[relative_summary$shift_type != "none", ]
  rows <- list(); idx <- 0L
  for (mode in unique(df$score_unit_mode)) {
    for (target in targets) {
      target_col <- paste0("obj_", target)
      for (bm in benchmarks) {
        if (target == bm) next
        bm_col <- paste0("obj_", bm)
        if (!(target_col %in% names(df)) || !(bm_col %in% names(df))) next
        for (st in unique(df$shift_type)) {
          tmp <- df[df$score_unit_mode == mode & df$shift_type == st, ]
          ok <- is.finite(tmp[[target_col]]) & is.finite(tmp[[bm_col]])
          idx <- idx + 1L
          rows[[idx]] <- data.frame(
            score_unit_mode = mode,
            target_policy = target,
            benchmark = bm,
            shift_type = st,
            n_cells = sum(ok),
            cell_wins = sum(tmp[[target_col]][ok] < tmp[[bm_col]][ok]),
            cell_win_rate = ifelse(sum(ok) > 0, sum(tmp[[target_col]][ok] < tmp[[bm_col]][ok]) / sum(ok), NA_real_),
            mean_relative = mean(tmp[[target_col]][ok] / tmp[[bm_col]][ok], na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, rows)
}

make_replicate_win_summary <- function(all_runs) {
  id_cols <- c("score_unit_mode", "kappa", "shift_prob", "shift_type", "sim_id")
  wide <- reshape(all_runs[, c(id_cols, "policy", "objective")],
                  idvar = id_cols, timevar = "policy", direction = "wide")
  names(wide) <- gsub("objective\\.", "obj_", names(wide))

  targets <- c("debt_threshold_tuned", "debt_utility_tuned", "hybrid_utility_tuned")
  benchmarks <- c("calendar10", "calendar_tuned", "cusum_tuned", "alarm_tuned",
                  "debt_threshold_tuned", "always_retrain", "never_retrain")
  rows <- list(); idx <- 0L
  for (mode in unique(wide$score_unit_mode)) {
    for (target in targets) {
      target_col <- paste0("obj_", target)
      for (bm in benchmarks) {
        if (target == bm) next
        bm_col <- paste0("obj_", bm)
        if (!(target_col %in% names(wide)) || !(bm_col %in% names(wide))) next
        for (k in unique(wide$kappa)) {
          for (st in unique(wide$shift_type)) {
            tmp <- wide[wide$score_unit_mode == mode & wide$kappa == k & wide$shift_type == st, ]
            ok <- is.finite(tmp[[target_col]]) & is.finite(tmp[[bm_col]])
            idx <- idx + 1L
            rows[[idx]] <- data.frame(
              score_unit_mode = mode,
              kappa = k,
              shift_type = st,
              target_policy = target,
              benchmark = bm,
              n_rep = sum(ok),
              replicate_win_rate = ifelse(sum(ok) > 0, mean(tmp[[target_col]][ok] < tmp[[bm_col]][ok]), NA_real_),
              mean_target = mean(tmp[[target_col]][ok], na.rm = TRUE),
              mean_benchmark = mean(tmp[[bm_col]][ok], na.rm = TRUE),
              relative_mean = mean(tmp[[target_col]][ok], na.rm = TRUE) / mean(tmp[[bm_col]][ok], na.rm = TRUE),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  do.call(rbind, rows)
}

make_score_unit_sensitivity_summary <- function(cell_win_summary) {
  # Compact table for the manuscript/supplement: how robust are the headline
  # wins to the median, q75, and mean score-unit definitions?
  df <- cell_win_summary[cell_win_summary$shift_type != "none", ]
  if (nrow(df) == 0L) return(data.frame())
  key <- paste(df$score_unit_mode, df$target_policy, df$benchmark, sep = "|")
  groups <- split(df, key)
  rows <- lapply(groups, function(g) {
    data.frame(
      score_unit_mode = g$score_unit_mode[1],
      target_policy = g$target_policy[1],
      benchmark = g$benchmark[1],
      n_cells = sum(g$n_cells, na.rm = TRUE),
      cell_wins = sum(g$cell_wins, na.rm = TRUE),
      cell_win_rate = safe_ratio(sum(g$cell_wins, na.rm = TRUE), sum(g$n_cells, na.rm = TRUE)),
      mean_relative = mean(g$mean_relative, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}



make_primary_table1_summary <- function(all_runs, score_unit_mode = "q75", n_boot = 5000L,
                                        seed = 20260429L) {
  # Headline table used in the manuscript. The relative objective is computed
  # at the scenario-cell level: first average each policy over the matched
  # held-out paths in a cell, then take the target/benchmark ratio, then average
  # those ratios over the 72 non-stable cells. Bootstrap intervals are paired
  # Monte Carlo intervals conditional on the synthetic v13d design.
  df <- all_runs[all_runs$score_unit_mode == score_unit_mode &
                   all_runs$shift_type != "none", ]
  if (nrow(df) == 0L) return(data.frame())

  id_cols <- c("kappa", "shift_prob", "shift_type", "sim_id")
  wide <- reshape(df[, c(id_cols, "policy", "objective")],
                  idvar = id_cols, timevar = "policy", direction = "wide")
  names(wide) <- gsub("^objective\\.", "", names(wide))
  cell_key <- paste(wide$kappa, wide$shift_prob, wide$shift_type, sep = "|")
  cells <- unique(cell_key)

  pairs <- data.frame(
    target_label = c("Age-adjusted debt threshold", "Age-adjusted debt threshold",
                     "Debt utility", "Debt utility",
                     "Hybrid utility", "Hybrid utility"),
    target_policy = c("debt_threshold_tuned", "debt_threshold_tuned",
                      "debt_utility_tuned", "debt_utility_tuned",
                      "hybrid_utility_tuned", "hybrid_utility_tuned"),
    benchmark_label = c("Tuned calendar", "Tuned CUSUM",
                        "Tuned calendar", "Tuned CUSUM",
                        "Tuned calendar", "Tuned CUSUM"),
    benchmark_policy = c("calendar_tuned", "cusum_tuned",
                         "calendar_tuned", "cusum_tuned",
                         "calendar_tuned", "cusum_tuned"),
    stringsAsFactors = FALSE
  )

  cell_ratio <- function(target, benchmark) {
    out <- numeric(length(cells))
    for (i in seq_along(cells)) {
      g <- wide[cell_key == cells[i], ]
      out[i] <- mean(g[[target]], na.rm = TRUE) / mean(g[[benchmark]], na.rm = TRUE)
    }
    out
  }

  set.seed(seed)
  boot_ci <- function(target, benchmark) {
    vals <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
      s <- 0
      for (cell in cells) {
        g <- wide[cell_key == cell, ]
        idx <- sample(seq_len(nrow(g)), size = nrow(g), replace = TRUE)
        s <- s + mean(g[[target]][idx], na.rm = TRUE) / mean(g[[benchmark]][idx], na.rm = TRUE)
      }
      vals[b] <- s / length(cells)
    }
    stats::quantile(vals, c(0.025, 0.975), names = FALSE, type = 7)
  }

  rows <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    target <- pairs$target_policy[i]
    benchmark <- pairs$benchmark_policy[i]
    rel <- cell_ratio(target, benchmark)
    ci <- boot_ci(target, benchmark)

    rows[[i]] <- data.frame(
      score_unit_mode = score_unit_mode,
      target_label = pairs$target_label[i],
      target_policy = target,
      benchmark_label = pairs$benchmark_label[i],
      benchmark_policy = benchmark,
      n_cells = length(rel),
      cell_wins = sum(rel < 1),
      wins_label = paste0(sum(rel < 1), "/", length(rel)),
      mean_relative = mean(rel),
      ci_025 = ci[1],
      ci_975 = ci[2],
      median_relative = median(rel),
      iqr_25 = stats::quantile(rel, 0.25, names = FALSE, type = 7),
      iqr_75 = stats::quantile(rel, 0.75, names = FALSE, type = 7),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

make_summaries <- function(all_runs) {
  policy_summary <- summarise_by_group(
    all_runs,
    c("score_unit_mode", "score_unit", "kappa", "lambda", "shift_prob", "shift_type", "policy")
  )
  relative_summary <- make_relative_summary(policy_summary)
  cell_win_summary <- make_cell_win_summary(relative_summary)
  replicate_win_summary <- make_replicate_win_summary(all_runs)
  score_unit_sensitivity <- make_score_unit_sensitivity_summary(cell_win_summary)
  primary_table1_q75 <- make_primary_table1_summary(all_runs, score_unit_mode = "q75")
  list(policy_summary = policy_summary,
       relative_summary = relative_summary,
       cell_win_summary = cell_win_summary,
       replicate_win_summary = replicate_win_summary,
       score_unit_sensitivity = score_unit_sensitivity,
       primary_table1_q75 = primary_table1_q75)
}

# =============================================================================
# 9. Figures, optional ggplot2
# =============================================================================

make_figures <- function(policy_summary, relative_summary, cal_df_adj, illustrative_series) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    log_msg("ggplot2 not installed. Skipping figures.")
    return(invisible(NULL))
  }
  library(ggplot2)

  rel_primary <- relative_summary[relative_summary$score_unit_mode == "q75", ]
  rel_long <- data.frame()
  add_rel <- function(df, name, label) {
    if (!(name %in% names(df))) return(data.frame())
    data.frame(
      score_unit_mode = df$score_unit_mode,
      kappa = df$kappa,
      shift_prob = df$shift_prob,
      shift_type = df$shift_type,
      comparison = label,
      relative_objective = df[[name]],
      stringsAsFactors = FALSE
    )
  }
  rel_long <- rbind(
    add_rel(rel_primary, "rel_debt_utility_vs_calendar", "Debt utility vs tuned calendar"),
    add_rel(rel_primary, "rel_hybrid_utility_vs_calendar", "Hybrid utility vs tuned calendar"),
    add_rel(rel_primary, "rel_debt_utility_vs_cusum", "Debt utility vs tuned CUSUM"),
    add_rel(rel_primary, "rel_hybrid_utility_vs_cusum", "Hybrid utility vs tuned CUSUM")
  )
  rel_long <- rel_long[rel_long$shift_type != "none" & is.finite(rel_long$relative_objective), ]
  rel_long$relative_plot <- pmax(rel_long$relative_objective, 1e-3)

  if (nrow(rel_long) > 0L) {
    p1 <- ggplot(rel_long[rel_long$comparison %in% c("Debt utility vs tuned calendar", "Hybrid utility vs tuned calendar"), ],
                 aes(x = factor(kappa), y = relative_plot, fill = factor(shift_prob))) +
      geom_col(position = "dodge", width = 0.7, color = "white") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      facet_grid(comparison ~ shift_type) +
      scale_y_log10() +
      labs(title = "Utility-calibrated debt policies relative to tuned calendar",
           subtitle = "75th-percentile score-unit scaling. Values below one favor the debt policy.",
           x = "Cost ratio kappa", y = "Relative objective, log scale", fill = "Shift probability") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")
    ggsave(file.path(OUTPUT_DIR, "fig1_relative_to_tuned_calendar_v13d.pdf"), p1, width = 12, height = 7)
  }

  p2_df <- policy_summary[policy_summary$score_unit_mode == "q75" & policy_summary$kappa == 1.0 & policy_summary$shift_type != "none", ]
  if (nrow(p2_df) > 0L) {
    p2_mean <- aggregate(mean_objective ~ shift_type + policy, data = p2_df, FUN = mean)
    p2 <- ggplot(p2_mean, aes(x = reorder(policy, mean_objective), y = mean_objective)) +
      geom_col(width = 0.7) + facet_wrap(~shift_type, scales = "free_y") + coord_flip() +
      labs(title = "Policy comparison at kappa = 1",
           subtitle = "75th-percentile score-unit scaling, averaged over shift probabilities.",
           x = NULL, y = "Mean predictive-regret objective") +
      theme_minimal(base_size = 11)
    ggsave(file.path(OUTPUT_DIR, "fig2_policy_comparison_kappa1_v13d.pdf"), p2, width = 10, height = 6)
  }

  p3_df <- policy_summary[policy_summary$score_unit_mode == "q75" & policy_summary$shift_type == "none", ]
  if (nrow(p3_df) > 0L) {
    p3_mean <- aggregate(mean_retrains ~ kappa + policy, data = p3_df, FUN = mean)
    p3 <- ggplot(p3_mean, aes(x = factor(kappa), y = mean_retrains, fill = policy)) +
      geom_col(position = "dodge", width = 0.7) +
      labs(title = "No-shift retraining behavior",
           subtitle = "75th-percentile score-unit scaling. Always and never baselines are included.",
           x = "Cost ratio kappa", y = "Mean retrains") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")
    ggsave(file.path(OUTPUT_DIR, "fig3_no_shift_retrains_v13d.pdf"), p3, width = 11, height = 6)
  }

  z <- cal_df_adj$debt_z
  br <- unique(stats::quantile(z[is.finite(z)], probs = seq(0, 1, length.out = 21), na.rm = TRUE))
  if (length(br) > 2L) {
    bins <- cut(z, breaks = br, include.lowest = TRUE)
    bdf <- aggregate(cbind(debt_z, regret_eval) ~ bins, data = cal_df_adj, FUN = mean)
    p4 <- ggplot(bdf, aes(x = debt_z, y = regret_eval)) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      labs(title = "Calibration: age-adjusted debt versus predictive regret",
           x = "Age-adjusted debt z-score", y = "Mean one-period predictive regret") +
      theme_minimal(base_size = 11)
    ggsave(file.path(OUTPUT_DIR, "fig4_calibration_debt_regret_v13d.pdf"), p4, width = 8, height = 5)
  }

  if (!is.null(illustrative_series) && nrow(illustrative_series) > 0L) {
    p5 <- ggplot(illustrative_series, aes(x = t)) +
      geom_line(aes(y = regret_eval), linewidth = 0.7) +
      geom_line(aes(y = rhat), linewidth = 0.7, linetype = "dashed") +
      geom_hline(aes(yintercept = lambda), linetype = "dotted") +
      geom_vline(data = illustrative_series[illustrative_series$action == 1L, ], aes(xintercept = t), alpha = 0.5) +
      labs(title = "Illustrative hybrid-utility decision path",
           subtitle = "Solid = realized evaluation regret, dashed = calibrated expected regret, dotted = retraining cost.",
           x = "Monitoring period", y = "Predictive-regret units") +
      theme_minimal(base_size = 11)
    ggsave(file.path(OUTPUT_DIR, "fig5_illustrative_decision_path_v13d.pdf"), p5, width = 10, height = 5.5)
  }

  sens_rows <- data.frame()
  add_sens <- function(name, label) {
    if (!(name %in% names(relative_summary))) return(data.frame())
    tmp <- relative_summary[relative_summary$shift_type != "none", ]
    data.frame(score_unit_mode = tmp$score_unit_mode,
               shift_type = tmp$shift_type,
               comparison = label,
               relative_objective = tmp[[name]],
               stringsAsFactors = FALSE)
  }
  sens_rows <- rbind(
    add_sens("rel_debt_threshold_vs_calendar", "Debt threshold vs calendar"),
    add_sens("rel_debt_threshold_vs_cusum", "Debt threshold vs CUSUM"),
    add_sens("rel_hybrid_utility_vs_calendar", "Hybrid vs calendar"),
    add_sens("rel_hybrid_utility_vs_cusum", "Hybrid vs CUSUM"),
    add_sens("rel_debt_utility_vs_calendar", "Debt utility vs calendar"),
    add_sens("rel_debt_utility_vs_cusum", "Debt utility vs CUSUM")
  )
  sens_rows <- sens_rows[is.finite(sens_rows$relative_objective), ]
  if (nrow(sens_rows) > 0L) {
    sens_mean <- aggregate(relative_objective ~ score_unit_mode + shift_type + comparison, data = sens_rows, FUN = mean)
    p6 <- ggplot(sens_mean, aes(x = score_unit_mode, y = relative_objective, fill = comparison)) +
      geom_col(position = "dodge", width = 0.7) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      facet_wrap(~shift_type) +
      labs(title = "Score-unit sensitivity",
           subtitle = "Mean relative objective across non-stable cells. Values below one favor the debt policy.",
           x = "Score-unit definition", y = "Mean relative objective") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")
    ggsave(file.path(OUTPUT_DIR, "fig6_score_unit_sensitivity_v13d.pdf"), p6, width = 12, height = 6)
  }
}

# =============================================================================
# 10. Main workflow
# =============================================================================

log_msg("=== Learning Debt v13d expanded-grid sensitivity simulation ===")
log_msg("Run mode: %s", RUN_MODE)
log_msg("Output directory: %s", OUTPUT_DIR)
log_msg("Resume: %s", ifelse(RESUME, "yes", "no"))
log_msg("T=%d, N_INIT=%d, update=%d, monitor=%d, eval=%d", T_PERIODS, N_INIT, N_UPDATE, N_MONITOR, N_EVAL)
log_msg("Paths per cell: calibration=%d, tuning=%d, test=%d", N_CAL_PATHS_PER_CELL, N_TUNE_PATHS_PER_CELL, N_TEST_PATHS_PER_CELL)

calibration_rows <- collect_calibration_rows(N_CAL_PATHS_PER_CELL)
age_adjuster <- fit_age_adjuster(calibration_rows)
calibration_rows_adj <- add_calibrated_features(calibration_rows, age_adjuster)
utility_models <- fit_utility_models(calibration_rows_adj)

score_units_df <- build_score_units(calibration_rows_adj)
candidate_grid <- build_candidate_grid(calibration_rows_adj)
tune_out <- tune_policies(candidate_grid, score_units_df, age_adjuster, utility_models, N_TUNE_PATHS_PER_CELL)
boundary_summary <- make_boundary_summary(tune_out$tuned_params, candidate_grid)
all_runs <- run_test_grid(tune_out$tuned_params, score_units_df, age_adjuster, utility_models, N_TEST_PATHS_PER_CELL)
summary_out <- make_summaries(all_runs)

# Illustrative path under the primary 75th-percentile score-unit mode.
set.seed(make_seed(SEED_ILLUS_BASE, "abrupt_coef", 0.05, 1L))
illustrative_path <- generate_path(T_PERIODS, 0.05, "abrupt_coef")
primary_mode <- "q75"
primary_score_unit <- score_units_df$score_unit[score_units_df$score_unit_mode == primary_mode][1]
if (!is.finite(primary_score_unit)) primary_score_unit <- score_units_df$score_unit[1]
illustrative_lambda <- 1.0 * primary_score_unit
illustrative_scale <- get_tuned_value(tune_out$tuned_params, primary_mode, 1.0, "hybrid_utility", 1)
illustrative_res <- simulate_policy(illustrative_path, "hybrid_utility", illustrative_lambda,
                                    age_adjuster, utility_models, illustrative_scale, TRUE)
illustrative_series <- illustrative_res$series
illustrative_series$shift_type <- "abrupt_coef"
illustrative_series$shift_prob <- 0.05
illustrative_series$kappa <- 1.0
illustrative_series$score_unit_mode <- primary_mode

# Save outputs.
write.csv(calibration_rows_adj, checkpoint_file("calibration_rows_v13d.csv"), row.names = FALSE)
write.csv(score_units_df, checkpoint_file("score_units_v13d.csv"), row.names = FALSE)
write.csv(candidate_grid, checkpoint_file("candidate_grid_v13d.csv"), row.names = FALSE)
write.csv(tune_out$tuning_summary, checkpoint_file("tuning_summary_v13d.csv"), row.names = FALSE)
write.csv(tune_out$tuned_params, checkpoint_file("tuned_parameters_v13d.csv"), row.names = FALSE)
write.csv(boundary_summary, checkpoint_file("boundary_summary_v13d.csv"), row.names = FALSE)
write.csv(summary_out$policy_summary, checkpoint_file("policy_summary_v13d.csv"), row.names = FALSE)
write.csv(summary_out$relative_summary, checkpoint_file("relative_summary_v13d.csv"), row.names = FALSE)
write.csv(summary_out$cell_win_summary, checkpoint_file("cell_win_summary_v13d.csv"), row.names = FALSE)
write.csv(summary_out$replicate_win_summary, checkpoint_file("replicate_win_summary_v13d.csv"), row.names = FALSE)
write.csv(summary_out$score_unit_sensitivity, checkpoint_file("score_unit_sensitivity_summary_v13d.csv"), row.names = FALSE)
write.csv(summary_out$primary_table1_q75, checkpoint_file("table1_primary_q75_v13d.csv"), row.names = FALSE)
write.csv(illustrative_series, checkpoint_file("illustrative_series_v13d.csv"), row.names = FALSE)

writeLines(capture.output(summary(age_adjuster$model)), checkpoint_file("age_adjuster_summary_v13d.txt"))
writeLines(capture.output(summary(utility_models$debt_model)), checkpoint_file("debt_utility_model_summary_v13d.txt"))
writeLines(capture.output(summary(utility_models$hybrid_model)), checkpoint_file("hybrid_utility_model_summary_v13d.txt"))
writeLines(capture.output(sessionInfo()), checkpoint_file("session_info_v13d.txt"))

make_figures(summary_out$policy_summary, summary_out$relative_summary, calibration_rows_adj, illustrative_series)

readme_lines <- c(
  "Learning Debt simulation outputs, v13d",
  "",
  paste0("Run mode: ", RUN_MODE),
  paste0("T_PERIODS: ", T_PERIODS),
  paste0("Calibration paths per scenario cell: ", N_CAL_PATHS_PER_CELL),
  paste0("Tuning paths per scenario cell: ", N_TUNE_PATHS_PER_CELL),
  paste0("Test paths per scenario cell: ", N_TEST_PATHS_PER_CELL),
  "",
  "Design summary:",
  "- Warm-started deployed and shadow posteriors.",
  "- Separate update, monitor, and evaluation batches.",
  "- Period-t actions affect period t+1 and later.",
  "- Objective is lambda times retrain count plus accumulated positive evaluation score gap.",
  "- Exact KL between NIG shadow and deployed posteriors is the debt signal.",
  "- Debt is age-adjusted using stable no-shift calibration paths.",
  "- Utility-calibrated debt and hybrid policies are tuned on calibration paths and evaluated on held-out paths.",
  "- Calendar, CUSUM, alarm, debt-threshold, always-retrain, and never-retrain baselines are included.",
  "- v13d expands the grids because v13c placed some baselines at boundary values.",
  "- Score-unit sensitivity is run for median, q75, and mean positive predictive-regret units.",
  "- The script checkpoints calibration, tuning, and held-out test paths and can resume after interruption.",
  "",
  "Expanded grids:",
  paste0("- Calendar intervals: ", paste(CAL_INTERVALS, collapse = ", ")),
  paste0("- CUSUM thresholds: ", paste(CUSUM_THRESHOLDS, collapse = ", ")),
  paste0("- Alarm quantiles: ", paste(ALARM_QUANTILES, collapse = ", ")),
  paste0("- Utility scales: ", paste(UTILITY_SCALES, collapse = ", ")),
  paste0("- Debt-threshold quantiles: ", paste(DEBT_THRESHOLD_PROBS, collapse = ", ")),
  "",
  "Main files:",
  "- all_runs_v13d.csv",
  "- policy_summary_v13d.csv",
  "- relative_summary_v13d.csv",
  "- cell_win_summary_v13d.csv",
  "- replicate_win_summary_v13d.csv",
  "- score_unit_sensitivity_summary_v13d.csv",
  "- table1_primary_q75_v13d.csv",
  "- tuning_summary_v13d.csv",
  "- tuned_parameters_v13d.csv",
  "- boundary_summary_v13d.csv",
  "- calibration_rows_v13d.csv",
  "- score_units_v13d.csv",
  "- session_info_v13d.txt"
)
writeLines(readme_lines, checkpoint_file("README_v13d.txt"))

log_msg("=== v13d complete. Outputs written to %s ===", OUTPUT_DIR)

cat("\n--- Score units ---\n")
print(score_units_df, row.names = FALSE)

cat("\n--- Tuned parameters ---\n")
print(tune_out$tuned_params[order(tune_out$tuned_params$score_unit_mode,
                                  tune_out$tuned_params$kappa,
                                  tune_out$tuned_params$policy), ], row.names = FALSE)

cat("\n--- Boundary summary ---\n")
print(boundary_summary[order(boundary_summary$score_unit_mode,
                             boundary_summary$kappa,
                             boundary_summary$policy), ], row.names = FALSE)

cat("\n--- Primary Table 1 summary, q75 score-unit mode ---\n")
print(summary_out$primary_table1_q75, row.names = FALSE)

cat("\n--- Cell win summary, non-stable cells ---\n")
print(summary_out$cell_win_summary[order(summary_out$cell_win_summary$score_unit_mode,
                                         summary_out$cell_win_summary$target_policy,
                                         summary_out$cell_win_summary$benchmark,
                                         summary_out$cell_win_summary$shift_type), ], row.names = FALSE)

cat("\n--- Mean relative objectives by score-unit mode and shift type ---\n")
rel <- summary_out$relative_summary[summary_out$relative_summary$shift_type != "none", ]
for (mode in unique(rel$score_unit_mode)) {
  for (st in unique(rel$shift_type)) {
    tmp <- rel[rel$score_unit_mode == mode & rel$shift_type == st, ]
    cat(sprintf("%s / %s: debt/calendar=%.3f, hybrid/calendar=%.3f, debt/CUSUM=%.3f, hybrid/CUSUM=%.3f\n",
                mode, st,
                mean(tmp$rel_debt_utility_vs_calendar, na.rm = TRUE),
                mean(tmp$rel_hybrid_utility_vs_calendar, na.rm = TRUE),
                mean(tmp$rel_debt_utility_vs_cusum, na.rm = TRUE),
                mean(tmp$rel_hybrid_utility_vs_cusum, na.rm = TRUE)))
  }
}
