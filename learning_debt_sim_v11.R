# =============================================================================
# Simulation Study: Learning Debt and Cost-Sensitive Bayesian Retraining
# v11 -- cadence sensitivity, uncertainty summaries, and one-period-lag robustness
#
# What this version adds beyond v10
# ---------------------------------
#  1. Calendar-cadence sensitivity over a grid of fixed retraining intervals.
#  2. A best-fixed-cadence benchmark selected ex post within that grid.
#  3. Replication counts and Monte Carlo uncertainty summaries for key ratios.
#  4. A one-period-lag robustness experiment for the online decision loop.
#  5. Output tables/figures/CSV files for all of the above.
#
# Important modeling disclosures
# ------------------------------
#  * The KL quantity implemented below is the exact KL divergence between two
#    normal-inverse-gamma posteriors in the conjugate simulation.
#  * c_churn and c_wait are one-step EXCESS LOSSES relative to the statewise
#    best action, not literal end-to-end operational spend.
#  * HMM emissions are calibrated to regime shift under a frozen deployment,
#    not directly to actionable staleness under an online policy.
#  * The proxy filter is not separately calibrated on its own signal
#    distribution; it maps a weighted z-score combination into the debt-filter
#    emission space.
#  * Same-period monitoring and action are used in the main experiment. A
#    one-period-lag robustness experiment is included below.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

set.seed(2026)

# =============================================================================
# 0. Global parameters
# =============================================================================

N_SIM      <- 300
T_PERIODS  <- 200
BURN_IN    <- 20
N_UPDATE   <- 5
N_HOLDOUT  <- 5
N_OBS      <- N_UPDATE + N_HOLDOUT

# Cost ratio grid: kappa = c_churn / c_wait
KAPPAS      <- c(0.1, 0.25, 0.5, 1.0, 2.0, 4.0)
SHIFT_PROBS <- c(0.02, 0.05, 0.10, 0.20)
SHIFT_TYPES <- c("none", "abrupt_coef", "variance", "gradual")

# Fixed calendar baselines
DEFAULT_CAL_INTERVAL <- 10L
CAL_INTERVALS        <- c(5L, 10L, 20L, 40L)
RUN_CADENCE_SENSITIVITY <- TRUE
CADENCE_N_SIM <- N_SIM

# One-period-lag robustness
RUN_LAG_ROBUSTNESS <- TRUE
LAG_ROBUST_N_SIM   <- 300
LAG_ROBUST_KAPPAS  <- c(0.5, 1.0, 2.0)
LAG_ROBUST_SHIFT_PROBS <- c(0.05, 0.10)
LAG_ROBUST_SHIFT_TYPES <- c("abrupt_coef", "gradual", "variance")

# Existing robustness experiments
RUN_ROBUSTNESS <- TRUE
ROBUST_N_SIM   <- 500
ROBUST_KAPPAS  <- c(0.5, 1.0, 2.0)
ROBUST_SHIFT_PROBS <- c(0.05, 0.10)
ROBUST_SHIFT_TYPES <- c("abrupt_coef", "gradual", "variance")
ROBUST_FIXED_HAZARD <- 0.05

# Filter backoff
HMM_P10_NO_SHIFT <- 0.00
HMM_P10_ABRUPT   <- 0.00
HMM_P10_VARIANCE <- 0.00
HMM_P10_GRADUAL  <- 0.02

# Materiality threshold for gradual drift staleness
STALE_THRESH <- 0.5

# =============================================================================
# 1. NIG prior and conjugate update
# =============================================================================

PRIOR <- list(mu_0 = 0, kappa_0 = 1, alpha_0 = 2, beta_0 = 1)

bayes_update <- function(prior, data) {
  if (is.null(data) || nrow(data) == 0) return(prior)

  ss_xy   <- sum(data$x * data$y)
  ss_xx   <- sum(data$x^2)
  n       <- nrow(data)
  kappa_n <- prior$kappa_0 + ss_xx
  mu_n    <- (prior$kappa_0 * prior$mu_0 + ss_xy) / kappa_n
  alpha_n <- prior$alpha_0 + n / 2
  beta_n  <- prior$beta_0 +
    0.5 * (sum(data$y^2) - kappa_n * mu_n^2 + prior$kappa_0 * prior$mu_0^2)

  list(mu_0 = mu_n, kappa_0 = kappa_n, alpha_0 = alpha_n, beta_0 = beta_n)
}

# =============================================================================
# 2. Exact KL learning debt in the conjugate NIG model
# =============================================================================

kl_nig_exact <- function(p, q) {
  ap <- p$alpha_0; bp <- p$beta_0; kp <- p$kappa_0; mp <- p$mu_0
  aq <- q$alpha_0; bq <- q$beta_0; kq <- q$kappa_0; mq <- q$mu_0

  if (ap <= 0 || bp <= 0 || kp <= 0 || aq <= 0 || bq <= 0 || kq <= 0) {
    return(NA_real_)
  }

  kl_ig <- aq * log(bp / bq) + lgamma(aq) - lgamma(ap) +
    (ap - aq) * digamma(ap) + ap * (bq / bp - 1)

  e_inv_sigma2_p <- ap / bp
  kl_beta_given_sigma2 <- 0.5 * (
    log(kp / kq) + kq / kp - 1 + kq * (mp - mq)^2 * e_inv_sigma2_p
  )

  max(kl_ig + kl_beta_given_sigma2, 0)
}

# =============================================================================
# 3. Holdout predictive log score
# =============================================================================

pred_log_score <- function(posterior, data) {
  if (is.null(data) || nrow(data) == 0) return(NA_real_)

  pred_mean <- posterior$mu_0 * data$x
  pred_var  <- (posterior$beta_0 / posterior$alpha_0) *
    (1 + data$x^2 / posterior$kappa_0)
  df <- 2 * posterior$alpha_0

  ls <- dt((data$y - pred_mean) / sqrt(pred_var), df = df, log = TRUE) -
    log(sqrt(pred_var))
  mean(ls)
}

# =============================================================================
# 4. Metrics and helpers
# =============================================================================

brier_score <- function(rho, Z) mean((rho - Z)^2)

brier_skill <- function(bs, bs_naive) {
  if (!is.finite(bs_naive) || bs_naive <= 1e-12) return(NA_real_)
  1 - bs / bs_naive
}

running_z <- function(x, s, ss, n) {
  if (n < 2) return(0)
  mn <- s / n
  vr <- max(ss / n - mn^2, 1e-8)
  (x - mn) / sqrt(vr)
}

safe_ratio <- function(num, den) ifelse(is.finite(den) & den > 0, num / den, NA_real_)

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(stats::quantile(x, prob = prob, names = FALSE, type = 7))
}

se_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

get_hmm_p10 <- function(shift_type) {
  switch(
    shift_type,
    "none"        = HMM_P10_NO_SHIFT,
    "abrupt_coef" = HMM_P10_ABRUPT,
    "variance"    = HMM_P10_VARIANCE,
    "gradual"     = HMM_P10_GRADUAL,
    HMM_P10_ABRUPT
  )
}

# =============================================================================
# 5. DGP: generate world regimes, parameters, and data path
# =============================================================================

generate_path <- function(T, shift_prob, shift_type) {
  G_true  <- integer(T)
  params  <- vector("list", T)
  cur_par <- list(beta = 0, sigma2 = 1)

  for (t in seq_len(T)) {
    prev_G <- if (t == 1) 0L else G_true[t - 1]

    if (shift_type == "none") {
      G_true[t] <- 0L

    } else if (shift_type == "abrupt_coef") {
      if (t > BURN_IN && prev_G == 0 && runif(1) < shift_prob) {
        G_true[t]    <- 1L
        cur_par$beta <- cur_par$beta + rnorm(1, 0, 2)
      } else {
        G_true[t] <- prev_G
      }

    } else if (shift_type == "variance") {
      if (t > BURN_IN && prev_G == 0 && runif(1) < shift_prob) {
        G_true[t]      <- 1L
        cur_par$sigma2 <- cur_par$sigma2 * runif(1, 3, 6)
      } else {
        G_true[t] <- prev_G
      }

    } else if (shift_type == "gradual") {
      if (t > BURN_IN && prev_G == 0 && runif(1) < shift_prob) {
        G_true[t] <- 1L
      } else {
        G_true[t] <- prev_G
      }
      if (G_true[t] == 1) {
        cur_par$beta <- cur_par$beta + rnorm(1, 0, 0.15)
      }
    }

    params[[t]] <- cur_par
  }

  obs <- lapply(seq_len(T), function(t) {
    x <- rnorm(N_OBS)
    y <- params[[t]]$beta * x + rnorm(N_OBS, sd = sqrt(params[[t]]$sigma2))
    data.frame(x = x, y = y)
  })

  list(G_true = G_true, params = params, obs = obs)
}

# =============================================================================
# 6. Policy-specific staleness definition
# =============================================================================

compute_staleness <- function(t, shift_type, G_true, last_retrain_t,
                              dep_post, shadow_post) {
  if (shift_type == "none") return(0L)

  if (shift_type %in% c("abrupt_coef", "variance")) {
    shift_onset <- which(G_true[seq_len(t)] == 1)[1]
    if (is.na(shift_onset)) return(0L)
    return(as.integer(shift_onset > last_retrain_t))
  }

  if (shift_type == "gradual") {
    prior_sd <- sqrt(PRIOR$beta_0 / PRIOR$alpha_0 / PRIOR$kappa_0)
    return(as.integer(abs(shadow_post$mu_0 - dep_post$mu_0) >
                        STALE_THRESH * prior_sd))
  }

  0L
}

# =============================================================================
# 7. One-step HMM update
# =============================================================================

hmm_update <- function(rho_prev, signal, shift_prob, p10,
                       em_mu_0, em_sd_0, em_mu_1, em_sd_1) {
  rho_pred <- shift_prob * (1 - rho_prev) + (1 - p10) * rho_prev
  lik_0    <- dnorm(signal, em_mu_0, em_sd_0)
  lik_1    <- dnorm(signal, em_mu_1, em_sd_1)
  denom    <- lik_0 * (1 - rho_pred) + lik_1 * rho_pred
  lik_1 * rho_pred / max(denom, 1e-12)
}

# =============================================================================
# 8. HMM emission calibration
# =============================================================================

calibrate_hmm_emissions <- function(shift_type, n_cal = 100,
                                    T = T_PERIODS, shift_prob = 0.05) {
  if (shift_type == "none") {
    cat(sprintf("  %s: neutral emissions\n", shift_type))
    return(list(mu_0 = 0.5, sd_0 = 0.5, mu_1 = 0.5, sd_1 = 0.5))
  }

  cat(sprintf("  Calibrating emissions for %s...\n", shift_type))
  d_stable  <- numeric(0)
  d_shifted <- numeric(0)

  for (i in seq_len(n_cal)) {
    path <- generate_path(T, shift_prob, shift_type)
    dep_post <- PRIOR
    shadow_post <- PRIOR

    for (t in seq_len(T)) {
      obs_u <- path$obs[[t]][seq_len(N_UPDATE), , drop = FALSE]
      shadow_post <- bayes_update(shadow_post, obs_u)
      d_val <- kl_nig_exact(shadow_post, dep_post)

      if (path$G_true[t] == 0) {
        d_stable  <- c(d_stable, d_val)
      } else {
        d_shifted <- c(d_shifted, d_val)
      }
    }
  }

  ds <- sqrt(pmax(d_stable, 0))
  dh <- sqrt(pmax(d_shifted, 0))

  list(
    mu_0 = mean(ds),
    sd_0 = max(sd(ds), 0.01),
    mu_1 = if (length(dh) > 1) mean(dh) else mean(ds) + 1.0,
    sd_1 = if (length(dh) > 1) max(sd(dh), 0.01) else 0.5
  )
}

build_emissions <- function(mode = c("family", "pooled"), n_cal = 100,
                            T = T_PERIODS, shift_prob = 0.05) {
  mode <- match.arg(mode)

  if (mode == "family") {
    cat("Calibrating family-specific emissions...\n")
    out <- list()
    for (st in SHIFT_TYPES) {
      out[[st]] <- calibrate_hmm_emissions(st, n_cal = n_cal, T = T,
                                           shift_prob = shift_prob)
      cat(sprintf("    %s: mu_0=%.3f sd_0=%.3f mu_1=%.3f sd_1=%.3f\n",
                  st,
                  out[[st]]$mu_0, out[[st]]$sd_0,
                  out[[st]]$mu_1, out[[st]]$sd_1))
    }
    return(out)
  }

  cat("Calibrating pooled emissions across non-none shift families...\n")
  pooled_stable  <- numeric(0)
  pooled_shifted <- numeric(0)

  for (st in setdiff(SHIFT_TYPES, "none")) {
    cat(sprintf("  contributing %s...\n", st))
    for (i in seq_len(n_cal)) {
      path <- generate_path(T, shift_prob, st)
      dep_post <- PRIOR
      shadow_post <- PRIOR

      for (t in seq_len(T)) {
        obs_u <- path$obs[[t]][seq_len(N_UPDATE), , drop = FALSE]
        shadow_post <- bayes_update(shadow_post, obs_u)
        d_val <- kl_nig_exact(shadow_post, dep_post)

        if (path$G_true[t] == 0) {
          pooled_stable  <- c(pooled_stable, d_val)
        } else {
          pooled_shifted <- c(pooled_shifted, d_val)
        }
      }
    }
  }

  ds <- sqrt(pmax(pooled_stable, 0))
  dh <- sqrt(pmax(pooled_shifted, 0))
  pooled <- list(
    mu_0 = mean(ds),
    sd_0 = max(sd(ds), 0.01),
    mu_1 = mean(dh),
    sd_1 = max(sd(dh), 0.01)
  )

  out <- list()
  out[["none"]] <- list(mu_0 = 0.5, sd_0 = 0.5, mu_1 = 0.5, sd_1 = 0.5)
  for (st in setdiff(SHIFT_TYPES, "none")) out[[st]] <- pooled
  cat(sprintf("    pooled: mu_0=%.3f sd_0=%.3f mu_1=%.3f sd_1=%.3f\n",
              pooled$mu_0, pooled$sd_0, pooled$mu_1, pooled$sd_1))
  out
}

EMISSIONS_FAMILY <- build_emissions("family", n_cal = 100)
EMISSIONS_POOLED <- build_emissions("pooled", n_cal = 50)

# =============================================================================
# 9. simulate_policy: run one policy on one pre-generated path
#
# decision_lag = 0: same-period monitoring and action (main experiment)
# decision_lag = 1: action chosen at t is applied at t+1 (lag robustness)
# =============================================================================

simulate_policy <- function(path, policy, kappa, shift_prob,
                            emissions = EMISSIONS_FAMILY,
                            hazard_override = NULL,
                            cal_interval = DEFAULT_CAL_INTERVAL,
                            cusum_thresh = 4,
                            decision_lag = 0L) {
  T          <- length(path$G_true)
  G_true     <- path$G_true
  shift_type <- attr(path, "shift_type")
  em         <- emissions[[shift_type]]
  p10        <- get_hmm_p10(shift_type)

  shift_prob_eff <- if (!is.null(hazard_override)) {
    hazard_override
  } else if (shift_type == "none") {
    0
  } else {
    shift_prob
  }

  c_churn   <- kappa / (1 + kappa)
  c_wait    <- 1 / (1 + kappa)
  threshold <- c_churn / (c_churn + c_wait)

  dep_post     <- PRIOR
  shadow_post  <- PRIOR
  last_retrain <- 0L

  rho_prev  <- 0
  cusum_val <- 0
  d_history <- numeric(0)

  sg_sum <- 0; sg_ss <- 0; sg_n <- 0
  pd_sum <- 0; pd_ss <- 0; pd_n <- 0
  re_sum <- 0; re_ss <- 0; re_n <- 0

  total_cost <- 0
  retrain    <- integer(T)      # actual executed retrains
  decisions  <- integer(T)      # actions chosen at time t
  Z_true     <- integer(T)
  rho_series <- rep(NA_real_, T)

  D_exact      <- numeric(T)
  score_gap    <- numeric(T)
  param_div    <- numeric(T)
  resid_exceed <- numeric(T)

  pending_action <- 0L
  pending_shadow_post <- NULL

  reset_deployment_relative_state <- function() {
    if (policy %in% c("debt_filter", "proxy_filter")) {
      rho_prev <<- 0
    }
    if (policy == "proxy_filter") {
      sg_sum <<- 0; sg_ss <<- 0; sg_n <<- 0
      pd_sum <<- 0; pd_ss <<- 0; pd_n <<- 0
      re_sum <<- 0; re_ss <<- 0; re_n <<- 0
    }
    if (policy == "cusum") {
      cusum_val <<- 0
    }
    if (policy == "alarm") {
      d_history <<- numeric(0)
    }
  }

  for (t in seq_len(T)) {
    executed_action <- 0L

    # Apply lagged retrain chosen last period, if any.
    if (decision_lag == 1L && pending_action == 1L) {
      dep_post <- pending_shadow_post
      last_retrain <- t - 1L
      executed_action <- 1L
      pending_action <- 0L
      pending_shadow_post <- NULL
      reset_deployment_relative_state()
    }

    obs_u <- path$obs[[t]][seq_len(N_UPDATE), , drop = FALSE]
    obs_h <- path$obs[[t]][(N_UPDATE + 1):N_OBS, , drop = FALSE]

    shadow_post <- bayes_update(shadow_post, obs_u)

    Z_t <- compute_staleness(t, shift_type, G_true, last_retrain,
                             dep_post, shadow_post)
    Z_true[t] <- Z_t

    sg <- pred_log_score(shadow_post, obs_h) - pred_log_score(dep_post, obs_h)
    if (is.na(sg)) sg <- 0
    pd <- abs(shadow_post$mu_0 - dep_post$mu_0)
    re <- mean(abs(obs_h$y - dep_post$mu_0 * obs_h$x) >
                 2 * sqrt(dep_post$beta_0 / dep_post$alpha_0))
    dp <- kl_nig_exact(shadow_post, dep_post)

    score_gap[t]    <- sg
    param_div[t]    <- pd
    resid_exceed[t] <- re
    D_exact[t]      <- dp

    sg_z <- running_z(sg, sg_sum, sg_ss, sg_n)
    pd_z <- running_z(pd, pd_sum, pd_ss, pd_n)
    re_z <- running_z(re, re_sum, re_ss, re_n)

    sg_sum <- sg_sum + sg; sg_ss <- sg_ss + sg^2; sg_n <- sg_n + 1
    pd_sum <- pd_sum + pd; pd_ss <- pd_ss + pd^2; pd_n <- pd_n + 1
    re_sum <- re_sum + re; re_ss <- re_ss + re^2; re_n <- re_n + 1

    alarm_thresh <- if (length(d_history) >= 10) safe_quantile(d_history, 0.90) else Inf
    d_history <- c(d_history, dp)

    cusum_val <- max(0, cusum_val + max(sg, 0))

    action_now <- switch(
      policy,
      "debt_filter" = {
        sqrt_dp <- sqrt(max(dp, 0))
        rho_new <- hmm_update(rho_prev, sqrt_dp, shift_prob_eff, p10,
                              em$mu_0, em$sd_0, em$mu_1, em$sd_1)
        rho_prev <- rho_new
        rho_series[t] <- rho_new
        as.integer(rho_new > threshold)
      },
      "proxy_filter" = {
        proxy_raw <- 0.5 * sg_z + 0.4 * pd_z + 0.2 * re_z
        proxy_scaled <- em$mu_0 + proxy_raw * em$sd_0
        rho_new <- hmm_update(rho_prev, proxy_scaled, shift_prob_eff, p10,
                              em$mu_0, em$sd_0, em$mu_1, em$sd_1)
        rho_prev <- rho_new
        rho_series[t] <- rho_new
        as.integer(rho_new > threshold)
      },
      "calendar" = as.integer(t %% cal_interval == 0),
      "cusum"    = as.integer(cusum_val > cusum_thresh),
      "alarm"    = as.integer(dp > alarm_thresh),
      stop("Unknown policy: ", policy)
    )
    decisions[t] <- action_now

    effective_action <- if (decision_lag == 0L) action_now else executed_action

    if (effective_action == 1L && Z_t == 0L) total_cost <- total_cost + c_churn
    if (effective_action == 0L && Z_t == 1L) total_cost <- total_cost + c_wait

    retrain[t] <- effective_action

    if (decision_lag == 0L) {
      if (action_now == 1L) {
        dep_post <- shadow_post
        last_retrain <- t
        reset_deployment_relative_state()
      }
    } else {
      if (action_now == 1L) {
        pending_action <- 1L
        pending_shadow_post <- shadow_post
      }
    }
  }

  list(
    cost            = total_cost,
    retrain         = retrain,
    decisions       = decisions,
    Z_true          = Z_true,
    rho_series      = rho_series,
    D_exact         = D_exact,
    score_gap       = score_gap,
    param_div       = param_div,
    resid_exceed    = resid_exceed,
    shift_prob_eff  = shift_prob_eff
  )
}

# =============================================================================
# 10. Single simulation replicate: main policy comparison on one path
# =============================================================================

run_one <- function(kappa, shift_prob, shift_type,
                    emissions = EMISSIONS_FAMILY,
                    hazard_override = NULL,
                    decision_lag = 0L,
                    cal_interval = DEFAULT_CAL_INTERVAL) {
  path <- generate_path(T_PERIODS, shift_prob, shift_type)
  attr(path, "shift_type") <- shift_type

  r_debt  <- simulate_policy(path, "debt_filter",  kappa, shift_prob,
                             emissions = emissions, hazard_override = hazard_override,
                             decision_lag = decision_lag, cal_interval = cal_interval)
  r_proxy <- simulate_policy(path, "proxy_filter", kappa, shift_prob,
                             emissions = emissions, hazard_override = hazard_override,
                             decision_lag = decision_lag, cal_interval = cal_interval)
  r_cal   <- simulate_policy(path, "calendar",     kappa, shift_prob,
                             emissions = emissions, hazard_override = hazard_override,
                             decision_lag = decision_lag, cal_interval = cal_interval)
  r_cusum <- simulate_policy(path, "cusum",        kappa, shift_prob,
                             emissions = emissions, hazard_override = hazard_override,
                             decision_lag = decision_lag, cal_interval = cal_interval)
  r_alarm <- simulate_policy(path, "alarm",        kappa, shift_prob,
                             emissions = emissions, hazard_override = hazard_override,
                             decision_lag = decision_lag, cal_interval = cal_interval)

  bs_naive_debt  <- brier_score(rep(r_debt$shift_prob_eff,  T_PERIODS), r_debt$Z_true)
  bs_naive_proxy <- brier_score(rep(r_proxy$shift_prob_eff, T_PERIODS), r_proxy$Z_true)
  bs_debt        <- brier_score(r_debt$rho_series,  r_debt$Z_true)
  bs_proxy       <- brier_score(r_proxy$rho_series, r_proxy$Z_true)

  bss_debt  <- brier_skill(bs_debt,  bs_naive_debt)
  bss_proxy <- brier_skill(bs_proxy, bs_naive_proxy)

  cor_sg <- if (sd(r_debt$D_exact) > 0 && sd(r_debt$score_gap) > 0) {
    cor(r_debt$D_exact, r_debt$score_gap, method = "spearman")
  } else {
    NA_real_
  }
  cor_pd <- if (sd(r_debt$D_exact) > 0 && sd(r_debt$param_div) > 0) {
    cor(r_debt$D_exact, r_debt$param_div, method = "spearman")
  } else {
    NA_real_
  }

  world_shift_onset <- which(path$G_true == 1)[1]
  policy_onset <- function(result) {
    if (shift_type %in% c("abrupt_coef", "variance")) {
      return(world_shift_onset)
    }
    which(result$Z_true == 1)[1]
  }
  delay_fn <- function(retrain_vec, onset) {
    if (is.na(onset)) return(NA_real_)
    first <- which(retrain_vec == 1 & seq_along(retrain_vec) >= onset)[1]
    if (is.na(first)) return(T_PERIODS - onset)
    first - onset
  }

  data.frame(
    kappa          = kappa,
    shift_prob     = shift_prob,
    shift_type     = shift_type,
    decision_lag   = decision_lag,
    cal_interval   = cal_interval,
    cost_debt      = r_debt$cost,
    cost_proxy     = r_proxy$cost,
    cost_calendar  = r_cal$cost,
    cost_cusum     = r_cusum$cost,
    cost_alarm     = r_alarm$cost,
    rel_debt_vs_calendar_rep = safe_ratio(r_debt$cost, r_cal$cost),
    rel_debt_vs_cusum_rep    = safe_ratio(r_debt$cost, r_cusum$cost),
    rel_debt_vs_alarm_rep    = safe_ratio(r_debt$cost, r_alarm$cost),
    bss_debt       = bss_debt,
    bss_proxy      = bss_proxy,
    cor_sg         = cor_sg,
    cor_pd         = cor_pd,
    delay_debt     = delay_fn(r_debt$retrain,  policy_onset(r_debt)),
    delay_calendar = delay_fn(r_cal$retrain,   policy_onset(r_cal)),
    delay_cusum    = delay_fn(r_cusum$retrain, policy_onset(r_cusum)),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 11. Calendar-cadence sensitivity on one path
# =============================================================================

run_one_cadence <- function(kappa, shift_prob, shift_type,
                            emissions = EMISSIONS_FAMILY,
                            hazard_override = NULL,
                            decision_lag = 0L,
                            cal_intervals = CAL_INTERVALS) {
  path <- generate_path(T_PERIODS, shift_prob, shift_type)
  attr(path, "shift_type") <- shift_type

  r_debt <- simulate_policy(path, "debt_filter", kappa, shift_prob,
                            emissions = emissions, hazard_override = hazard_override,
                            decision_lag = decision_lag, cal_interval = DEFAULT_CAL_INTERVAL)

  world_shift_onset <- which(path$G_true == 1)[1]
  debt_onset <- if (shift_type %in% c("abrupt_coef", "variance")) {
    world_shift_onset
  } else {
    which(r_debt$Z_true == 1)[1]
  }
  delay_fn <- function(retrain_vec, onset) {
    if (is.na(onset)) return(NA_real_)
    first <- which(retrain_vec == 1 & seq_along(retrain_vec) >= onset)[1]
    if (is.na(first)) return(T_PERIODS - onset)
    first - onset
  }

  bind_rows(lapply(cal_intervals, function(ci) {
    r_cal <- simulate_policy(path, "calendar", kappa, shift_prob,
                             emissions = emissions, hazard_override = hazard_override,
                             decision_lag = decision_lag, cal_interval = ci)
    cal_onset <- if (shift_type %in% c("abrupt_coef", "variance")) {
      world_shift_onset
    } else {
      which(r_cal$Z_true == 1)[1]
    }
    data.frame(
      kappa = kappa,
      shift_prob = shift_prob,
      shift_type = shift_type,
      decision_lag = decision_lag,
      calendar_interval = ci,
      cost_debt = r_debt$cost,
      cost_calendar = r_cal$cost,
      rel_debt_vs_calendar_rep = safe_ratio(r_debt$cost, r_cal$cost),
      delay_debt = delay_fn(r_debt$retrain, debt_onset),
      delay_calendar = delay_fn(r_cal$retrain, cal_onset),
      stringsAsFactors = FALSE
    )
  }))
}

# =============================================================================
# 12. Grid runners
# =============================================================================

run_grid <- function(grid, runner = c("main", "cadence"),
                     emissions = EMISSIONS_FAMILY,
                     hazard_override = NULL,
                     decision_lag = 0L,
                     progress_every = 2000) {
  runner <- match.arg(runner)
  out <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    if (i %% progress_every == 0) cat(sprintf("  %d / %d\n", i, nrow(grid)))
    r <- grid[i, ]
    out[[i]] <- tryCatch(
      if (runner == "main") {
        run_one(r$kappa, r$shift_prob, r$shift_type,
                emissions = emissions,
                hazard_override = hazard_override,
                decision_lag = decision_lag)
      } else {
        run_one_cadence(r$kappa, r$shift_prob, r$shift_type,
                        emissions = emissions,
                        hazard_override = hazard_override,
                        decision_lag = decision_lag)
      },
      error = function(e) {
        warning(sprintf("Run failed for kappa=%s, shift_prob=%s, shift_type=%s: %s",
                        r$kappa, r$shift_prob, r$shift_type, e$message))
        if (runner == "main") {
          data.frame(
            kappa = r$kappa,
            shift_prob = r$shift_prob,
            shift_type = r$shift_type,
            decision_lag = decision_lag,
            cal_interval = DEFAULT_CAL_INTERVAL,
            cost_debt = NA_real_,
            cost_proxy = NA_real_,
            cost_calendar = NA_real_,
            cost_cusum = NA_real_,
            cost_alarm = NA_real_,
            rel_debt_vs_calendar_rep = NA_real_,
            rel_debt_vs_cusum_rep = NA_real_,
            rel_debt_vs_alarm_rep = NA_real_,
            bss_debt = NA_real_,
            bss_proxy = NA_real_,
            cor_sg = NA_real_,
            cor_pd = NA_real_,
            delay_debt = NA_real_,
            delay_calendar = NA_real_,
            delay_cusum = NA_real_,
            stringsAsFactors = FALSE
          )
        } else {
          data.frame(
            kappa = r$kappa,
            shift_prob = r$shift_prob,
            shift_type = r$shift_type,
            decision_lag = decision_lag,
            calendar_interval = CAL_INTERVALS,
            cost_debt = NA_real_,
            cost_calendar = NA_real_,
            rel_debt_vs_calendar_rep = NA_real_,
            delay_debt = NA_real_,
            delay_calendar = NA_real_,
            stringsAsFactors = FALSE
          )
        }
      }
    )
  }
  bind_rows(out)
}

# =============================================================================
# 13. Main experiment
# =============================================================================

main_grid <- expand.grid(
  kappa      = KAPPAS,
  shift_prob = SHIFT_PROBS,
  shift_type = SHIFT_TYPES,
  sim_id     = seq_len(N_SIM),
  stringsAsFactors = FALSE
)

cat(sprintf("\nRunning main experiment: %d runs\n", nrow(main_grid)))
all_df <- run_grid(main_grid, runner = "main", emissions = EMISSIONS_FAMILY,
                   hazard_override = NULL, decision_lag = 0L)
cat("Main experiment complete.\n")

# =============================================================================
# 14. Main summaries and uncertainty summaries
# =============================================================================

main_summary <- all_df %>%
  group_by(kappa, shift_prob, shift_type) %>%
  summarise(
    n_rep          = sum(is.finite(cost_debt)),
    mean_debt      = mean(cost_debt,     na.rm = TRUE),
    mean_proxy     = mean(cost_proxy,    na.rm = TRUE),
    mean_calendar  = mean(cost_calendar, na.rm = TRUE),
    mean_cusum     = mean(cost_cusum,    na.rm = TRUE),
    mean_alarm     = mean(cost_alarm,    na.rm = TRUE),
    se_debt        = se_mean(cost_debt),
    se_calendar    = se_mean(cost_calendar),
    mean_bss_debt  = mean(bss_debt,      na.rm = TRUE),
    mean_bss_proxy = mean(bss_proxy,     na.rm = TRUE),
    mean_cor_sg    = mean(cor_sg,        na.rm = TRUE),
    mean_cor_pd    = mean(cor_pd,        na.rm = TRUE),
    mean_delay_debt     = mean(delay_debt,     na.rm = TRUE),
    mean_delay_calendar = mean(delay_calendar, na.rm = TRUE),
    mean_delay_cusum    = mean(delay_cusum,    na.rm = TRUE),
    pct_debt_best = mean(
      cost_debt <= pmin(cost_calendar, cost_cusum, cost_alarm, na.rm = TRUE),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    rel_debt_vs_calendar  = mean_debt  / mean_calendar,
    rel_debt_vs_cusum     = mean_debt  / mean_cusum,
    rel_debt_vs_alarm     = mean_debt  / mean_alarm,
    rel_proxy_vs_calendar = mean_proxy / mean_calendar
  )

uncertainty_summary <- all_df %>%
  group_by(kappa, shift_prob, shift_type) %>%
  summarise(
    n_ratio_calendar = sum(is.finite(rel_debt_vs_calendar_rep)),
    ratio_cal_mean   = mean(rel_debt_vs_calendar_rep, na.rm = TRUE),
    ratio_cal_se     = se_mean(rel_debt_vs_calendar_rep),
    ratio_cal_q025   = safe_quantile(rel_debt_vs_calendar_rep, 0.025),
    ratio_cal_median = safe_quantile(rel_debt_vs_calendar_rep, 0.500),
    ratio_cal_q975   = safe_quantile(rel_debt_vs_calendar_rep, 0.975),
    ratio_cusum_mean = mean(rel_debt_vs_cusum_rep, na.rm = TRUE),
    ratio_cusum_se   = se_mean(rel_debt_vs_cusum_rep),
    ratio_cusum_q025 = safe_quantile(rel_debt_vs_cusum_rep, 0.025),
    ratio_cusum_q975 = safe_quantile(rel_debt_vs_cusum_rep, 0.975),
    ratio_alarm_mean = mean(rel_debt_vs_alarm_rep, na.rm = TRUE),
    ratio_alarm_se   = se_mean(rel_debt_vs_alarm_rep),
    ratio_alarm_q025 = safe_quantile(rel_debt_vs_alarm_rep, 0.025),
    ratio_alarm_q975 = safe_quantile(rel_debt_vs_alarm_rep, 0.975),
    .groups = "drop"
  )

# =============================================================================
# 15. Calendar-cadence sensitivity
# =============================================================================

cadence_summary <- NULL
best_cadence_summary <- NULL
if (RUN_CADENCE_SENSITIVITY) {
  cadence_grid <- expand.grid(
    kappa      = KAPPAS,
    shift_prob = SHIFT_PROBS,
    shift_type = SHIFT_TYPES,
    sim_id     = seq_len(CADENCE_N_SIM),
    stringsAsFactors = FALSE
  )

  cat(sprintf("\nRunning calendar-cadence sensitivity: %d paths\n", nrow(cadence_grid)))
  cadence_df <- run_grid(cadence_grid, runner = "cadence", emissions = EMISSIONS_FAMILY,
                         hazard_override = NULL, decision_lag = 0L, progress_every = 1000)

  cadence_summary <- cadence_df %>%
    group_by(kappa, shift_prob, shift_type, calendar_interval) %>%
    summarise(
      n_rep = sum(is.finite(cost_debt)),
      mean_debt = mean(cost_debt, na.rm = TRUE),
      mean_calendar = mean(cost_calendar, na.rm = TRUE),
      se_calendar = se_mean(cost_calendar),
      rel_debt_vs_calendar = mean_debt / mean_calendar,
      ratio_rep_mean = mean(rel_debt_vs_calendar_rep, na.rm = TRUE),
      ratio_rep_se   = se_mean(rel_debt_vs_calendar_rep),
      ratio_rep_q025 = safe_quantile(rel_debt_vs_calendar_rep, 0.025),
      ratio_rep_q975 = safe_quantile(rel_debt_vs_calendar_rep, 0.975),
      mean_delay_debt = mean(delay_debt, na.rm = TRUE),
      mean_delay_calendar = mean(delay_calendar, na.rm = TRUE),
      .groups = "drop"
    )

  best_cadence_summary <- cadence_summary %>%
    group_by(kappa, shift_prob, shift_type) %>%
    summarise(
      mean_debt = first(mean_debt),
      best_calendar_interval = calendar_interval[which.min(mean_calendar)],
      mean_best_calendar = min(mean_calendar, na.rm = TRUE),
      rel_debt_vs_best_calendar = mean_debt / mean_best_calendar,
      .groups = "drop"
    )
}

# =============================================================================
# 16. Existing robustness experiments
# =============================================================================

robustness_summary <- NULL
if (RUN_ROBUSTNESS) {
  robust_grid <- expand.grid(
    kappa      = ROBUST_KAPPAS,
    shift_prob = ROBUST_SHIFT_PROBS,
    shift_type = ROBUST_SHIFT_TYPES,
    sim_id     = seq_len(ROBUST_N_SIM),
    stringsAsFactors = FALSE
  )

  cat(sprintf("\nRunning robustness experiment: %d runs per condition\n", nrow(robust_grid)))

  cat("  robustness A: default emissions + true hazard\n")
  rob_default <- run_grid(robust_grid, runner = "main", emissions = EMISSIONS_FAMILY,
                          hazard_override = NULL, decision_lag = 0L, progress_every = 1000) %>%
    mutate(robust_mode = "default")

  cat("  robustness B: pooled emissions + true hazard\n")
  rob_pooled <- run_grid(robust_grid, runner = "main", emissions = EMISSIONS_POOLED,
                         hazard_override = NULL, decision_lag = 0L, progress_every = 1000) %>%
    mutate(robust_mode = "pooled_emissions")

  cat("  robustness C: default emissions + fixed hazard\n")
  rob_hazard <- run_grid(robust_grid, runner = "main", emissions = EMISSIONS_FAMILY,
                         hazard_override = ROBUST_FIXED_HAZARD,
                         decision_lag = 0L, progress_every = 1000) %>%
    mutate(robust_mode = "fixed_hazard")

  robustness_df <- bind_rows(rob_default, rob_pooled, rob_hazard)

  robustness_summary <- robustness_df %>%
    group_by(robust_mode, kappa, shift_prob, shift_type) %>%
    summarise(
      n_rep = sum(is.finite(cost_debt)),
      mean_debt = mean(cost_debt, na.rm = TRUE),
      mean_calendar = mean(cost_calendar, na.rm = TRUE),
      mean_cusum = mean(cost_cusum, na.rm = TRUE),
      mean_alarm = mean(cost_alarm, na.rm = TRUE),
      mean_bss_debt = mean(bss_debt, na.rm = TRUE),
      rel_debt_vs_calendar = mean_debt / mean_calendar,
      rel_debt_vs_cusum    = mean_debt / mean_cusum,
      rel_debt_vs_alarm    = mean_debt / mean_alarm,
      .groups = "drop"
    )
}

# =============================================================================
# 17. One-period-lag robustness
# =============================================================================

lag_robustness_summary <- NULL
if (RUN_LAG_ROBUSTNESS) {
  lag_grid <- expand.grid(
    kappa      = LAG_ROBUST_KAPPAS,
    shift_prob = LAG_ROBUST_SHIFT_PROBS,
    shift_type = LAG_ROBUST_SHIFT_TYPES,
    sim_id     = seq_len(LAG_ROBUST_N_SIM),
    stringsAsFactors = FALSE
  )

  cat(sprintf("\nRunning one-period-lag robustness: %d runs\n", nrow(lag_grid)))
  lag_df <- run_grid(lag_grid, runner = "main", emissions = EMISSIONS_FAMILY,
                     hazard_override = NULL, decision_lag = 1L, progress_every = 1000)

  lag_robustness_summary <- lag_df %>%
    group_by(kappa, shift_prob, shift_type) %>%
    summarise(
      n_rep = sum(is.finite(cost_debt)),
      mean_debt = mean(cost_debt, na.rm = TRUE),
      mean_calendar = mean(cost_calendar, na.rm = TRUE),
      mean_cusum = mean(cost_cusum, na.rm = TRUE),
      mean_alarm = mean(cost_alarm, na.rm = TRUE),
      mean_bss_debt = mean(bss_debt, na.rm = TRUE),
      mean_delay_debt = mean(delay_debt, na.rm = TRUE),
      rel_debt_vs_calendar = mean_debt / mean_calendar,
      rel_debt_vs_cusum    = mean_debt / mean_cusum,
      rel_debt_vs_alarm    = mean_debt / mean_alarm,
      .groups = "drop"
    )

  lag_comparison_summary <- lag_robustness_summary %>%
    left_join(
      main_summary %>%
        filter(kappa %in% LAG_ROBUST_KAPPAS,
               shift_prob %in% LAG_ROBUST_SHIFT_PROBS,
               shift_type %in% LAG_ROBUST_SHIFT_TYPES) %>%
        select(kappa, shift_prob, shift_type,
               rel_debt_vs_calendar_nolag = rel_debt_vs_calendar,
               rel_debt_vs_cusum_nolag = rel_debt_vs_cusum,
               mean_bss_debt_nolag = mean_bss_debt),
      by = c("kappa", "shift_prob", "shift_type")
    ) %>%
    mutate(
      delta_rel_calendar = rel_debt_vs_calendar - rel_debt_vs_calendar_nolag,
      delta_rel_cusum    = rel_debt_vs_cusum - rel_debt_vs_cusum_nolag,
      delta_bss          = mean_bss_debt - mean_bss_debt_nolag
    )
} else {
  lag_comparison_summary <- NULL
}

# =============================================================================
# 18. Illustrative run for figures
# =============================================================================

set.seed(42)
illus_path <- generate_path(T_PERIODS, 0.05, "abrupt_coef")
attr(illus_path, "shift_type") <- "abrupt_coef"

illus_debt  <- simulate_policy(illus_path, "debt_filter",  1.0, 0.05)
illus_proxy <- simulate_policy(illus_path, "proxy_filter", 1.0, 0.05)
illus_cal   <- simulate_policy(illus_path, "calendar",     1.0, 0.05)

illus_df <- data.frame(
  t            = seq_len(T_PERIODS),
  G_true       = illus_path$G_true,
  Z_debt       = illus_debt$Z_true,
  rho_debt     = illus_debt$rho_series,
  rho_proxy    = illus_proxy$rho_series,
  D_exact      = illus_debt$D_exact,
  threshold    = 0.5,
  retrain_debt = illus_debt$retrain,
  retrain_cal  = illus_cal$retrain
)

shift_start <- min(which(illus_df$G_true == 1), na.rm = TRUE)
if (!is.finite(shift_start)) shift_start <- T_PERIODS

# =============================================================================
# 19. Figures
# =============================================================================

fig1_data <- main_summary %>%
  filter(shift_type != "none", is.finite(rel_debt_vs_calendar), rel_debt_vs_calendar > 0) %>%
  mutate(rel_plot = pmax(rel_debt_vs_calendar, 1e-3))

p1 <- ggplot(fig1_data,
             aes(x = factor(kappa), y = rel_plot, fill = factor(shift_prob))) +
  geom_col(position = "dodge", width = 0.7, color = "white") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  facet_wrap(~shift_type, ncol = 2,
             labeller = labeller(shift_type = c(
               abrupt_coef = "Abrupt coefficient",
               variance    = "Variance shift",
               gradual     = "Gradual drift"))) +
  scale_fill_brewer(palette = "Blues", name = "Shift probability") +
  scale_y_log10(labels = label_number(accuracy = 0.1)) +
  labs(
    title    = "Debt-filter excess loss relative to calendar retraining",
    subtitle = "Log scale used so variance-shift failures remain visible; values below 1 favor debt-filter",
    x = expression(paste("Cost ratio ", kappa, " = ", c[churn], "/", c[wait])),
    y = "Relative excess loss (debt-filter / calendar, log scale)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p2 <- main_summary %>%
  filter(shift_type == "abrupt_coef", kappa == 1.0) %>%
  select(shift_prob, rel_debt_vs_calendar, rel_proxy_vs_calendar,
         rel_debt_vs_cusum, rel_debt_vs_alarm) %>%
  pivot_longer(-shift_prob, names_to = "metric", values_to = "rel") %>%
  mutate(label = recode(metric,
    rel_debt_vs_calendar  = "Debt-filter vs Calendar",
    rel_proxy_vs_calendar = "Proxy-filter vs Calendar",
    rel_debt_vs_cusum     = "Debt-filter vs CUSUM",
    rel_debt_vs_alarm     = "Debt-filter vs Alarm-only"
  )) %>%
  ggplot(aes(x = factor(shift_prob), y = rel, color = label, group = label)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  scale_color_brewer(palette = "Set1", name = "") +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  labs(title = "Policy comparisons: abrupt shift, symmetric excess-loss weights",
       x = "Shift probability", y = "Relative excess loss ratio") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p3a <- illus_df %>%
  ggplot(aes(x = t)) +
  annotate("rect", xmin = shift_start, xmax = T_PERIODS, ymin = 0, ymax = 1,
           fill = "orange", alpha = 0.15) +
  geom_line(aes(y = rho_debt),  color = "steelblue", linewidth = 0.9) +
  geom_line(aes(y = rho_proxy), color = "tomato", linewidth = 0.7, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray30", linewidth = 0.8) +
  geom_vline(data = illus_df[illus_df$retrain_debt == 1, ],
             aes(xintercept = t), color = "steelblue", alpha = 0.5, linewidth = 0.4) +
  scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  labs(title = "Decision diagram: abrupt coefficient shift",
       subtitle = "Orange = shifted regime; blue solid = debt-filter; red dashed = proxy-filter",
       x = NULL, y = expression(hat(rho)[t])) +
  theme_minimal(base_size = 11)

p3b <- illus_df %>%
  ggplot(aes(x = t)) +
  annotate("rect", xmin = shift_start, xmax = T_PERIODS,
           ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.10) +
  geom_line(aes(y = as.integer(Z_debt)), color = "darkblue", linewidth = 0.8) +
  scale_y_continuous(limits = c(-0.1, 1.1)) +
  labs(subtitle = "Deployment staleness Z_t (debt-filter policy)",
       x = "Monitoring period", y = expression(Z[t])) +
  theme_minimal(base_size = 11)

p3 <- p3a / p3b + plot_layout(heights = c(2, 1))

p4 <- main_summary %>%
  filter(shift_type != "none") %>%
  select(kappa, shift_prob, shift_type, mean_cor_sg, mean_cor_pd) %>%
  pivot_longer(c(mean_cor_sg, mean_cor_pd), names_to = "stat", values_to = "cor") %>%
  mutate(stat_label = recode(stat,
    mean_cor_sg = "Score gap",
    mean_cor_pd = "Parameter divergence"
  )) %>%
  group_by(shift_prob, shift_type, stat_label) %>%
  summarise(cor = mean(cor, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = factor(shift_prob), y = cor, fill = stat_label)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~shift_type, ncol = 3) +
  scale_fill_manual(values = c("steelblue", "darkgreen"), name = "Monitoring statistic") +
  labs(title = "Spearman correlation: monitoring diagnostics vs exact-KL learning debt",
       subtitle = "Bars are averaged across cost-ratio settings kappa",
       x = "Shift probability", y = "Spearman correlation") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p5 <- main_summary %>%
  filter(shift_type == "abrupt_coef") %>%
  select(kappa, shift_prob, mean_delay_debt, mean_delay_calendar, mean_delay_cusum) %>%
  pivot_longer(c(mean_delay_debt, mean_delay_calendar, mean_delay_cusum),
               names_to = "policy", values_to = "delay") %>%
  mutate(policy_label = recode(policy,
    mean_delay_debt     = "Debt-filter",
    mean_delay_calendar = "Calendar",
    mean_delay_cusum    = "CUSUM"
  )) %>%
  ggplot(aes(x = factor(kappa), y = delay, fill = policy_label)) +
  geom_col(position = "dodge", width = 0.7, color = "white") +
  facet_wrap(~shift_prob, nrow = 1,
             labeller = labeller(shift_prob = function(x) paste0("p=", x))) +
  scale_fill_manual(values = c("steelblue", "gray50", "tomato"), name = "Policy") +
  labs(title = "Mean delay to retraining after first shift",
       subtitle = "Longer delay at high kappa reflects rational excess-loss weighting",
       x = expression(paste("Cost ratio ", kappa)),
       y = "Mean periods to first retrain") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p6 <- NULL
if (!is.null(robustness_summary)) {
  p6 <- robustness_summary %>%
    filter(shift_type == "abrupt_coef", kappa == 1.0) %>%
    ggplot(aes(x = factor(shift_prob), y = rel_debt_vs_calendar,
               color = robust_mode, group = robust_mode)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
    scale_color_brewer(palette = "Dark2", name = "Robustness mode") +
    scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    labs(title = "Robustness: debt-filter vs calendar in abrupt shifts",
         subtitle = "Compare default calibration, pooled emissions, and fixed hazard",
         x = "Shift probability", y = "Relative excess loss ratio") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

p7 <- NULL
if (!is.null(cadence_summary)) {
  p7 <- cadence_summary %>%
    filter(shift_type == "abrupt_coef", kappa == 1.0) %>%
    ggplot(aes(x = factor(calendar_interval), y = rel_debt_vs_calendar,
               color = factor(shift_prob), group = factor(shift_prob))) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
    scale_color_brewer(palette = "Set2", name = "Shift probability") +
    labs(title = "Calendar-cadence sensitivity: abrupt shifts, kappa = 1",
         subtitle = "Debt-filter relative to fixed cadences 5, 10, 20, and 40 periods",
         x = "Calendar interval", y = "Relative excess loss ratio") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

p8 <- NULL
if (!is.null(lag_comparison_summary)) {
  p8 <- lag_comparison_summary %>%
    filter(shift_type == "abrupt_coef", kappa == 1.0) %>%
    ggplot(aes(x = factor(shift_prob), y = delta_rel_calendar,
               fill = shift_type)) +
    geom_col(width = 0.7, color = "white") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
    scale_fill_brewer(palette = "Pastel1", guide = "none") +
    labs(title = "One-period-lag robustness: change in debt vs calendar ratio",
         subtitle = "Positive values mean lagged deployment makes the debt-filter less competitive",
         x = "Shift probability", y = expression(Delta~"relative excess loss")) +
    theme_minimal(base_size = 11)
}

# =============================================================================
# 20. Save outputs
# =============================================================================

output_dir <- "~/learning_debt_sim_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(output_dir, "fig1_cost_by_shift_type_v11.pdf"), p1, width = 10, height = 8, create.dir = TRUE)
ggsave(file.path(output_dir, "fig2_policy_comparisons_v11.pdf"), p2, width = 10, height = 5.5)
ggsave(file.path(output_dir, "fig3_decision_diagram_v11.pdf"), p3, width = 10, height = 8)
ggsave(file.path(output_dir, "fig4_dt_correlation_v11.pdf"), p4, width = 11, height = 5)
ggsave(file.path(output_dir, "fig5_delay_v11.pdf"), p5, width = 11, height = 5)
if (!is.null(p6)) ggsave(file.path(output_dir, "fig6_robustness_v11.pdf"), p6, width = 9, height = 5)
if (!is.null(p7)) ggsave(file.path(output_dir, "fig7_calendar_sensitivity_v11.pdf"), p7, width = 9, height = 5)
if (!is.null(p8)) ggsave(file.path(output_dir, "fig8_lag_robustness_v11.pdf"), p8, width = 8, height = 5)

write.csv(main_summary, file.path(output_dir, "main_summary_v11.csv"), row.names = FALSE)
write.csv(all_df, file.path(output_dir, "all_runs_v11.csv"), row.names = FALSE)
write.csv(uncertainty_summary, file.path(output_dir, "uncertainty_summary_v11.csv"), row.names = FALSE)
write.csv(illus_df, file.path(output_dir, "illustrative_run_v11.csv"), row.names = FALSE)
if (!is.null(robustness_summary)) {
  write.csv(robustness_summary, file.path(output_dir, "robustness_summary_v11.csv"), row.names = FALSE)
}
if (!is.null(cadence_summary)) {
  write.csv(cadence_summary, file.path(output_dir, "calendar_cadence_summary_v11.csv"), row.names = FALSE)
  write.csv(best_cadence_summary, file.path(output_dir, "best_calendar_cadence_summary_v11.csv"), row.names = FALSE)
}
if (!is.null(lag_robustness_summary)) {
  write.csv(lag_robustness_summary, file.path(output_dir, "lag_robustness_summary_v11.csv"), row.names = FALSE)
  write.csv(lag_comparison_summary, file.path(output_dir, "lag_comparison_summary_v11.csv"), row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info_v11.txt"))

readme_lines <- c(
  "Learning Debt simulation outputs (v11)",
  "",
  "Main additions relative to v10:",
  "- uncertainty_summary_v11.csv: Monte Carlo uncertainty for key relative-loss ratios",
  "- calendar_cadence_summary_v11.csv: sensitivity to fixed calendar cadences",
  "- best_calendar_cadence_summary_v11.csv: best fixed cadence within the tested grid",
  "- lag_robustness_summary_v11.csv: one-period-lag robustness summaries",
  "- lag_comparison_summary_v11.csv: change from no-lag to one-period-lag results",
  "",
  "Key modeling choices:",
  "- Exact KL divergence between NIG posteriors is used as the simulation debt quantity.",
  "- c_churn and c_wait are one-step excess losses relative to the statewise best action.",
  "- Same-period monitoring and action are used in the main experiment.",
  "- The proxy filter reuses emissions calibrated on the exact-KL debt signal.",
  "- Calendar sensitivity is evaluated over intervals 5, 10, 20, and 40 periods.",
  "- One-period-lag robustness is run on a reduced scenario grid.",
  "",
  "Main files:",
  "- main_summary_v11.csv: aggregated Monte Carlo summaries",
  "- all_runs_v11.csv: replicate-level results",
  "- uncertainty_summary_v11.csv: uncertainty summaries for ratios",
  "- calendar_cadence_summary_v11.csv: cadence sensitivity summaries",
  "- best_calendar_cadence_summary_v11.csv: best fixed cadence summaries",
  "- robustness_summary_v11.csv: pooled-emission and fixed-hazard robustness summaries",
  "- lag_robustness_summary_v11.csv: one-period-lag robustness summaries",
  "- session_info_v11.txt: R and package versions"
)
writeLines(readme_lines, file.path(output_dir, "README_v11.txt"))

cat("\n=== Done. Outputs written to", output_dir, "===\n")

# =============================================================================
# 21. Console tables
# =============================================================================

cat("\n--- No-shift sanity check: mean excess losses by policy ---\n")
main_summary %>%
  filter(shift_type == "none") %>%
  select(kappa, shift_prob, mean_debt, mean_proxy, mean_calendar, mean_cusum, mean_alarm) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  print(n = Inf)

cat("\n--- Debt-filter: relative excess loss vs benchmarks (abrupt shift) ---\n")
main_summary %>%
  filter(shift_type == "abrupt_coef") %>%
  select(kappa, shift_prob, rel_debt_vs_calendar, rel_debt_vs_cusum, rel_debt_vs_alarm,
         pct_debt_best, mean_delay_debt, mean_delay_calendar) %>%
  mutate(across(starts_with("rel"), ~round(.x, 3)),
         pct_debt_best = percent(pct_debt_best, accuracy = 1),
         across(starts_with("mean_delay"), ~round(.x, 1))) %>%
  print(n = Inf)

cat("\n--- Uncertainty summary: abrupt shift, debt vs calendar ---\n")
uncertainty_summary %>%
  filter(shift_type == "abrupt_coef") %>%
  select(kappa, shift_prob, ratio_cal_mean, ratio_cal_se, ratio_cal_q025, ratio_cal_q975) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  print(n = Inf)

cat("\n--- BSS and exact-KL learning-debt correlations (abrupt shift) ---\n")
main_summary %>%
  filter(shift_type == "abrupt_coef") %>%
  select(kappa, shift_prob, mean_bss_debt, mean_bss_proxy, mean_cor_sg, mean_cor_pd) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  print(n = Inf)

cat("\n--- Gradual drift: relative excess loss vs benchmarks ---\n")
main_summary %>%
  filter(shift_type == "gradual") %>%
  select(kappa, shift_prob, rel_debt_vs_calendar, rel_debt_vs_cusum, rel_debt_vs_alarm,
         mean_delay_debt) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  print(n = Inf)

cat("\n--- Variance shift: relative excess loss vs benchmarks ---\n")
main_summary %>%
  filter(shift_type == "variance") %>%
  select(kappa, shift_prob, rel_debt_vs_calendar, rel_debt_vs_cusum, rel_debt_vs_alarm,
         mean_delay_debt) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  print(n = Inf)

if (!is.null(best_cadence_summary)) {
  cat("\n--- Best fixed calendar cadence within {5,10,20,40} periods ---\n")
  best_cadence_summary %>%
    filter(shift_type == "abrupt_coef") %>%
    select(kappa, shift_prob, best_calendar_interval, rel_debt_vs_best_calendar) %>%
    mutate(across(where(is.numeric), ~round(.x, 3))) %>%
    print(n = Inf)
}

if (!is.null(robustness_summary)) {
  cat("\n--- Robustness: abrupt shift, kappa = 1 ---\n")
  robustness_summary %>%
    filter(shift_type == "abrupt_coef", kappa == 1.0) %>%
    select(robust_mode, shift_prob, rel_debt_vs_calendar, rel_debt_vs_cusum, mean_bss_debt) %>%
    mutate(across(where(is.numeric), ~round(.x, 3))) %>%
    print(n = Inf)
}

if (!is.null(lag_comparison_summary)) {
  cat("\n--- One-period-lag robustness: change in debt-filter relative losses ---\n")
  lag_comparison_summary %>%
    filter(shift_type == "abrupt_coef") %>%
    select(kappa, shift_prob, rel_debt_vs_calendar_nolag, rel_debt_vs_calendar,
           delta_rel_calendar, rel_debt_vs_cusum_nolag, rel_debt_vs_cusum, delta_rel_cusum) %>%
    mutate(across(where(is.numeric), ~round(.x, 3))) %>%
    print(n = Inf)
}
