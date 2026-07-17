#!/usr/bin/env Rscript
# =============================================================================
# verify_paper_numbers.R
# Learning Debt replication package: manuscript number verification
#
# Recomputes every number reported in the manuscript "Cost-sensitive
# retraining via posterior learning debt" from local replication outputs and
# prints PASS / FAIL / SKIP for each value. Base R only; no packages required.
#
# Folder resolution: for each dataset the script looks, in order, under
# replication/ in this repository, then under LD_OUTPUT_ROOT (default: the
# home directory). For the positive-accounting primary run it also falls back
# to the legacy outputs/ directory committed to this repository.
#
# Usage:
#   Rscript verify_paper_numbers.R
#   LD_OUTPUT_ROOT=/path/to/run/folders Rscript verify_paper_numbers.R
#
# Exit status is 0 when every available check passes and 1 otherwise, so the
# script can gate a release. Missing folders skip their group rather than
# fail, and the skip is reported.
# =============================================================================

options(stringsAsFactors = FALSE)

root_env <- path.expand(Sys.getenv("LD_OUTPUT_ROOT", unset = "~"))

find_dir <- function(leaf, extra = character(0)) {
  for (p in c(file.path("replication", leaf), file.path(root_env, leaf), extra)) {
    if (dir.exists(p)) return(p)
  }
  NA_character_
}

DIR_SIGNED   <- find_dir("learning_debt_sim_outputs_v13d_paper_signed")
DIR_POSITIVE <- find_dir("learning_debt_sim_outputs_v13d_paper_positive", extra = "outputs")
DIR_OMC_WELL <- find_dir("learning_debt_sim_outputs_v13d_paper_signed_outermc")
DIR_OMC_T4   <- find_dir("learning_debt_sim_outputs_v13d_paper_signed_student_t_outermc")
DIR_FRONTIER <- find_dir("learning_debt_sim_outputs_v13d_paper_signed_frontier")

cat("Folder resolution:\n")
show_dir <- function(nm, d) cat(sprintf("  %-34s %s\n", nm, ifelse(is.na(d), "NOT FOUND", d)))
show_dir("signed primary (grid):", DIR_SIGNED)
show_dir("positive primary (grid):", DIR_POSITIVE)
show_dir("selection audit, well specified:", DIR_OMC_WELL)
show_dir("selection audit, Student-t df 4:", DIR_OMC_T4)
show_dir("frontier sweep:", DIR_FRONTIER)
cat("\n")

read_if <- function(dir, leaf) {
  if (is.na(dir)) return(NULL)
  f <- file.path(dir, leaf)
  if (!file.exists(f)) return(NULL)
  read.csv(f)
}

RES <- list()
add_res <- function(group, name, status, detail = "") {
  RES[[length(RES) + 1L]] <<- data.frame(group = group, name = name,
                                         status = status, detail = detail)
}
chk <- function(group, name, computed, expected, tol = 6e-4) {
  if (is.null(computed) || length(computed) == 0L || is.na(computed[1])) {
    add_res(group, name, "FAIL", "computed value missing")
    return(invisible())
  }
  ok <- abs(as.numeric(computed[1]) - as.numeric(expected)) <= tol + 1e-12
  add_res(group, name, if (ok) "PASS" else "FAIL",
          sprintf("computed %.6g, expected %.6g", as.numeric(computed[1]), as.numeric(expected)))
}
skip_group <- function(group, n_checks, why) {
  add_res(group, sprintf("(%d checks)", n_checks), "SKIP", why)
}

# =============================================================================
# A. Score units (signed primary run; Section 4)
# =============================================================================
grp <- "A score units"
su <- read_if(DIR_SIGNED, "score_units_v13d.csv")
if (is.null(su)) {
  skip_group(grp, 3L, "score_units_v13d.csv not found in the signed primary folder")
} else {
  unit <- function(m) su$score_unit[su$score_unit_mode == m][1]
  chk(grp, "median unit", unit("median"), 0.0125, 6e-5)
  chk(grp, "q75 unit",    unit("q75"),    0.0807, 6e-5)
  chk(grp, "mean unit",   unit("mean"),   0.1490, 6e-5)
}

# =============================================================================
# B. Table 1: primary comparison, signed accounting, q75 unit (Section 5.1)
# =============================================================================
grp <- "B Table 1 (signed, q75)"
t1 <- read_if(DIR_SIGNED, "table1_primary_q75_v13d.csv")
T1_EXPECTED <- list(
  list("debt_threshold_tuned", "calendar_tuned", 72, 0.663, 0.660, 0.666, 0.660, 0.631, 0.692),
  list("debt_threshold_tuned", "cusum_tuned",    61, 0.967, 0.965, 0.969, 0.957, 0.941, 0.977),
  list("debt_utility_tuned",   "calendar_tuned", 66, 0.781, 0.775, 0.789, 0.774, 0.689, 0.857),
  list("debt_utility_tuned",   "cusum_tuned",     6, 1.143, 1.135, 1.154, 1.074, 1.035, 1.204),
  list("hybrid_utility_tuned", "calendar_tuned", 71, 0.780, 0.775, 0.787, 0.800, 0.688, 0.867),
  list("hybrid_utility_tuned", "cusum_tuned",     2, 1.141, 1.134, 1.150, 1.089, 1.054, 1.202)
)
if (is.null(t1)) {
  skip_group(grp, 42L, "table1_primary_q75_v13d.csv not found in the signed primary folder")
} else {
  for (e in T1_EXPECTED) {
    r <- t1[t1$target_policy == e[[1]] & t1$benchmark_policy == e[[2]], ]
    tag <- paste(sub("_tuned$", "", e[[1]]), "vs", sub("_tuned$", "", e[[2]]))
    if (nrow(r) != 1L) { add_res(grp, tag, "FAIL", "row not found in table1 CSV"); next }
    chk(grp, paste(tag, "wins"),   r$cell_wins,       e[[3]], 0)
    chk(grp, paste(tag, "mean"),   r$mean_relative,   e[[4]])
    chk(grp, paste(tag, "ci lo"),  r$ci_025,          e[[5]])
    chk(grp, paste(tag, "ci hi"),  r$ci_975,          e[[6]])
    chk(grp, paste(tag, "median"), r$median_relative, e[[7]])
    chk(grp, paste(tag, "iqr 25"), r$iqr_25,          e[[8]])
    chk(grp, paste(tag, "iqr 75"), r$iqr_75,          e[[9]])
  }
}

# =============================================================================
# Shared helper: debt threshold vs a benchmark from a relative summary
# =============================================================================
rel_stats <- function(rel, mode, bench_col) {
  s <- rel[rel$shift_type != "none" & rel$score_unit_mode == mode, ]
  t <- s$obj_debt_threshold_tuned
  b <- s[[bench_col]]
  ok <- is.finite(t) & is.finite(b) & b > 0
  list(wins = sum(t[ok] < b[ok]), ratio = mean(t[ok] / b[ok]), n = sum(ok))
}

# =============================================================================
# C. Cost-scale sensitivity, signed accounting (Section 5.2)
# =============================================================================
grp <- "C sensitivity (signed)"
rels <- read_if(DIR_SIGNED, "relative_summary_v13d.csv")
if (is.null(rels)) {
  skip_group(grp, 8L, "relative_summary_v13d.csv not found in the signed primary folder")
} else {
  for (e in list(list("median", 28, 1.037), list("q75", 61, 0.967), list("mean", 69, 0.958))) {
    st <- rel_stats(rels, e[[1]], "obj_cusum_tuned")
    chk(grp, paste("vs CUSUM", e[[1]], "wins"),  st$wins,  e[[2]], 0)
    chk(grp, paste("vs CUSUM", e[[1]], "ratio"), st$ratio, e[[3]])
  }
  st <- rel_stats(rels, "q75", "obj_calendar_tuned")
  chk(grp, "vs calendar q75 wins",  st$wins,  72, 0)
  chk(grp, "vs calendar q75 ratio", st$ratio, 0.663)
}

# =============================================================================
# D. Cost-scale sensitivity, positive-part accounting (Section 5.2)
# =============================================================================
grp <- "D sensitivity (positive part)"
relp <- read_if(DIR_POSITIVE, "relative_summary_v13d.csv")
if (is.null(relp)) {
  skip_group(grp, 12L, "relative_summary_v13d.csv not found in the positive folder or outputs/")
} else {
  for (e in list(list("median", "obj_cusum_tuned",    26, 1.049),
                 list("q75",    "obj_cusum_tuned",    58, 0.975),
                 list("mean",   "obj_cusum_tuned",    63, 0.965),
                 list("median", "obj_calendar_tuned", 48, 0.895),
                 list("q75",    "obj_calendar_tuned", 72, 0.677),
                 list("mean",   "obj_calendar_tuned", 72, 0.639))) {
    st <- rel_stats(relp, e[[1]], e[[2]])
    lbl <- paste(sub("obj_", "vs ", sub("_tuned$", "", e[[2]])), e[[1]])
    chk(grp, paste(lbl, "wins"),  st$wins,  e[[3]], 0)
    chk(grp, paste(lbl, "ratio"), st$ratio, e[[4]])
  }
}

# =============================================================================
# E. Table 2: stable no-shift behavior, positive part, q75 unit (Section 5.3)
# =============================================================================
grp <- "E no-shift table"
psp <- read_if(DIR_POSITIVE, "policy_summary_v13d.csv")
if (is.null(psp)) {
  skip_group(grp, 10L, "policy_summary_v13d.csv not found in the positive folder or outputs/")
} else {
  ns <- function(pol) {
    s <- psp[psp$shift_type == "none" & psp$score_unit_mode == "q75" & psp$policy == pol, ]
    c(retrains = mean(s$mean_retrains), objective = mean(s$mean_objective))
  }
  for (e in list(list("debt_threshold_tuned",   1.39,  0.565),
                 list("cusum_tuned",            8.06,  0.550),
                 list("calendar_tuned",        46.83,  2.773),
                 list("never_retrain",          0.00,  0.880),
                 list("always_retrain",       199.00, 21.064))) {
    v <- ns(e[[1]])
    chk(grp, paste(e[[1]], "retrains"),  v["retrains"],  e[[2]], 6e-3)
    chk(grp, paste(e[[1]], "objective"), v["objective"], e[[3]])
  }
}

# =============================================================================
# F. Table 3: selection-uncertainty audit (Section 5.4)
# =============================================================================
grp <- "F selection audit"
OMC_EXPECTED <- list(
  wellspec  = list(dir = DIR_OMC_WELL,
                   median = c(46, 1.021, 1.002, 1.037),
                   q75    = c(86, 0.954, 0.945, 0.959),
                   mean   = c(96, 0.949, 0.944, 0.959)),
  student_t = list(dir = DIR_OMC_T4,
                   median = c(63, 0.980, 0.936, 1.032),
                   q75    = c(88, 0.921, 0.897, 0.962),
                   mean   = c(94, 0.912, 0.889, 0.956))
)
for (cell in names(OMC_EXPECTED)) {
  spec <- OMC_EXPECTED[[cell]]
  omc <- read_if(spec$dir, "outer_mc_summary_v13d.csv")
  if (is.null(omc)) {
    skip_group(grp, 12L, sprintf("outer_mc_summary_v13d.csv not found for the %s cell", cell))
    next
  }
  for (m in c("median", "q75", "mean")) {
    r <- omc[omc$score_unit_mode == m, ]
    e <- spec[[m]]
    if (nrow(r) != 1L) { add_res(grp, paste(cell, m), "FAIL", "row not found"); next }
    chk(grp, paste(cell, m, "win pct"),   100 * r$dt_cusum_win_mean, e[1], 0.51)
    chk(grp, paste(cell, m, "ratio"),     r$dt_cusum_ratio_mean,     e[2])
    chk(grp, paste(cell, m, "ratio q05"), r$dt_cusum_ratio_q05,      e[3])
    chk(grp, paste(cell, m, "ratio q95"), r$dt_cusum_ratio_q95,      e[4])
  }
}

# =============================================================================
# G. Table 4: misspecification severity sweep (Section 6.2)
# =============================================================================
grp <- "G frontier table"
fs <- read_if(DIR_FRONTIER, "frontier_summary_v13d.csv")
FRONTIER_EXPECTED <- list(
  list("tail",      Inf, 1.022, 0.955, 0.951, 0.246),
  list("tail",      8,   1.016, 0.945, 0.938, 0.260),
  list("tail",      6,   0.998, 0.932, 0.925, 0.267),
  list("tail",      5,   0.996, 0.934, 0.928, 0.271),
  list("tail",      4,   0.965, 0.911, 0.907, 0.275),
  list("tail",      3,   1.012, 0.963, 0.945, 0.277),
  list("hetero",    0,   1.022, 0.955, 0.950, 0.246),
  list("hetero",    0.5, 1.020, 0.947, 0.937, 0.266),
  list("hetero",    1,   1.010, 0.938, 0.927, 0.276),
  list("hetero",    1.5, 1.004, 0.931, 0.919, 0.283),
  list("hetero",    2,   0.997, 0.927, 0.914, 0.287),
  list("hetero",    4,   0.981, 0.917, 0.903, 0.296),
  list("nonlinear", 0,   1.018, 0.952, 0.950, 0.250),
  list("nonlinear", 0.2, 1.035, 0.961, 0.956, 0.257),
  list("nonlinear", 0.4, 1.165, 1.028, 1.008, 0.238),
  list("nonlinear", 0.8, 1.821, 1.442, 1.384, 0.192)
)
if (is.null(fs)) {
  skip_group(grp, 64L, "frontier_summary_v13d.csv not found in the frontier folder")
} else {
  frow <- function(ax, sev, m) fs[fs$axis == ax & fs$severity == sev & fs$score_unit_mode == m, ]
  modes3 <- c("median", "q75", "mean")
  for (e in FRONTIER_EXPECTED) {
    tag <- sprintf("%s sev %s", e[[1]], format(e[[2]]))
    for (i in seq_along(modes3)) {
      r <- frow(e[[1]], e[[2]], modes3[i])
      chk(grp, paste(tag, modes3[i], "ratio"),
          if (nrow(r) == 1L) r$dt_cusum_ratio_mean else NA_real_, e[[2 + i]])
    }
    rmed <- frow(e[[1]], e[[2]], "median")
    chk(grp, paste(tag, "diagnostic"),
        if (nrow(rmed) == 1L) rmed$frac_deployed_better else NA_real_, e[[6]])
  }
}

# =============================================================================
# H. Tail-axis median replication bands (Section 6.2, in text)
# =============================================================================
grp <- "H tail median bands"
TAIL_BANDS <- list(list(Inf, 1.007, 1.035), list(8, 1.000, 1.038), list(6, 0.979, 1.012),
                   list(5, 0.966, 1.036), list(4, 0.932, 0.991), list(3, 0.885, 1.362))
if (is.null(fs)) {
  skip_group(grp, 12L, "frontier_summary_v13d.csv not found in the frontier folder")
} else {
  for (e in TAIL_BANDS) {
    r <- fs[fs$axis == "tail" & fs$severity == e[[1]] & fs$score_unit_mode == "median", ]
    tag <- sprintf("nu=%s", format(e[[1]]))
    chk(grp, paste(tag, "q05"), if (nrow(r) == 1L) r$dt_cusum_ratio_q05 else NA_real_, e[[2]])
    chk(grp, paste(tag, "q95"), if (nrow(r) == 1L) r$dt_cusum_ratio_q95 else NA_real_, e[[3]])
  }
}

# =============================================================================
# I. Crossing statements (Section 6.2)
# =============================================================================
grp <- "I crossings"
cr <- read_if(DIR_FRONTIER, "frontier_crossings_v13d.csv")
if (is.null(cr) || is.null(fs)) {
  skip_group(grp, 8L, "frontier_crossings_v13d.csv or frontier_summary_v13d.csv not found")
} else {
  crow <- function(ax) cr[cr$axis == ax & cr$score_unit_mode == "median", ]
  h <- crow("hetero")
  chk(grp, "hetero median crosses",   as.numeric(isTRUE(h$crosses[1])), 1, 0)
  chk(grp, "hetero cross from 1.5",   h$cross_from_severity, 1.5, 1e-6)
  chk(grp, "hetero cross to 2",       h$cross_to_severity,   2,   1e-6)
  tl <- crow("tail")
  chk(grp, "tail median crosses",     as.numeric(isTRUE(tl$crosses[1])), 1, 0)
  chk(grp, "tail cross from nu=8",    tl$cross_from_severity, 8, 1e-6)
  chk(grp, "tail cross to nu=6",      tl$cross_to_severity,   6, 1e-6)
  nl <- crow("nonlinear")
  chk(grp, "curvature median never crosses", as.numeric(isTRUE(nl$crosses[1])), 0, 0)
  r4 <- fs[fs$axis == "tail" & fs$severity == 4 & fs$score_unit_mode == "median", ]
  chk(grp, "nu=4 median band below one",
      as.numeric(nrow(r4) == 1L && r4$dt_cusum_ratio_q95 < 1), 1, 0)
}

# =============================================================================
# Report
# =============================================================================
res <- do.call(rbind, RES)
n_pass <- sum(res$status == "PASS")
n_fail <- sum(res$status == "FAIL")
skips  <- res[res$status == "SKIP", , drop = FALSE]

cat("============================================================\n")
cat(sprintf("Verification: %d PASS, %d FAIL, %d group(s) skipped\n",
            n_pass, n_fail, nrow(skips)))
cat("============================================================\n\n")

if (nrow(skips) > 0L) {
  cat("Skipped groups (folders not found; run the corresponding study or\n")
  cat("point LD_OUTPUT_ROOT at the folder that contains it):\n")
  for (i in seq_len(nrow(skips))) {
    cat(sprintf("  [%s] %s: %s\n", skips$group[i], skips$name[i], skips$detail[i]))
  }
  cat("\n")
}

if (n_fail > 0L) {
  cat("Failures:\n")
  f <- res[res$status == "FAIL", , drop = FALSE]
  for (i in seq_len(nrow(f))) {
    cat(sprintf("  [%s] %s: %s\n", f$group[i], f$name[i], f$detail[i]))
  }
  cat("\n")
  quit(status = 1)
}

cat("All available checks pass.\n")
