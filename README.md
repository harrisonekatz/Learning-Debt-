# Learning Debt

Replication package for the manuscript **Cost-sensitive retraining via
posterior learning debt** (H. E. Katz). The study defines posterior learning
debt, the Kullback-Leibler divergence from a continuously updated shadow
posterior to the frozen deployed posterior, adjusted for deployment age, and
attaches a cost-sensitive retraining rule priced in predictive-score units.
The evidence has three parts: a deployment-calibrated primary experiment in a
conjugate model with exact posterior divergences, a selection-uncertainty
audit that repeats the entire calibrate, tune, and evaluate pipeline across
independent replications, and a misspecification severity sweep that maps the
regime in which the debt signal outperforms tuned change detection, together
with an observable diagnostic that flags that regime from the monitoring
stream.

## Contents

| Item | Purpose |
| --- | --- |
| `learning_debt_sim_sensitivity.R` | Single script that runs all three studies |
| `verify_paper_numbers.R` | Recomputes every number reported in the manuscript from local outputs (183 checks) and exits nonzero on any mismatch |
| `replication/` | Replicate-level and summary outputs accompanying the submission (see `replication/MANIFEST.md`) |
| `outputs/` | Original single-pass primary run under the positive-part regret accounting, kept for the record; the verification script reads it directly |
| `CITATION.cff` | Citation metadata |

## Requirements

Base R (4.0 or later). No packages are required to run the studies or the
verification. If `ggplot2` is installed, grid mode additionally writes the
appendix figures.

## One script, three modes

The script selects its mode from environment variables. All modes checkpoint
their long stages and resume after interruption.

**1. Frontier sweep (the default).**

```
Rscript learning_debt_sim_sensitivity.R
```

Sweeps three misspecification axes at increasing severity (Student-t tails,
heteroskedasticity, static curvature as the negative control), running the
full pipeline with independent replications at every point: eight per point on
the two drift axes, three to four on the control axis, with common random
numbers across severity levels. Writes `frontier_per_rep_v13d.csv`,
`frontier_summary_v13d.csv`, and `frontier_crossings_v13d.csv` under
`~/learning_debt_sim_outputs_v13d_paper_signed_frontier/`. Reproduces Table 4,
the tail-axis replication bands, the crossing statements, and the regime
diagnostic. Roughly one and a half to two days on a single core.

**2. Selection-uncertainty audit.**

```
LD_FRONTIER=0 Rscript learning_debt_sim_sensitivity.R
```

Repeats the entire calibrate, tune, and evaluate pipeline 15 times with
independent seed streams for the well-specified baseline and the Student-t
(four degrees of freedom) cell, keeping the tuning sample at its deliberately
small size. Every calibrated quantity, the cost units included, is
re-estimated within each replication. Writes `outer_mc_per_rep_v13d.csv` and
`outer_mc_summary_v13d.csv` under
`~/learning_debt_sim_outputs_v13d_paper_signed_outermc/` and
`~/learning_debt_sim_outputs_v13d_paper_signed_student_t_outermc/`.
Reproduces Table 3. Roughly a day.

**3. Single-pass grid.**

```
LD_FRONTIER=0 LD_OUTER_REPS=0 Rscript learning_debt_sim_sensitivity.R
```

Runs the {positive, signed} x {wellspec, student_t, heteroskedastic,
nonlinear} grid at full inner size, one output folder per cell. The signed
well-specified cell reproduces Table 1 and the signed sensitivity numbers; the
positive well-specified cell reproduces the positive-part robustness numbers
and Table 2, and matches the committed `outputs/` directory exactly, random
stream included. Roughly two hours per cell; restrict cells with
`LD_REGRET_MODE` and `LD_DGP_MODE`, for example:

```
LD_FRONTIER=0 LD_OUTER_REPS=0 LD_REGRET_MODE=signed LD_DGP_MODE=wellspec Rscript learning_debt_sim_sensitivity.R
```

**Wiring check.** `LD_RUN_MODE=smoke Rscript learning_debt_sim_sensitivity.R`
runs the active mode at toy sizes in a few minutes.

## Environment switches

| Variable | Default | Meaning |
| --- | --- | --- |
| `LD_RUN_MODE` | `paper` | `smoke`, `paper`, or `full` experiment sizes |
| `LD_FRONTIER` | `1` | `0` leaves frontier mode |
| `LD_FRONTIER_REPS` | `8` | Replications per severity point |
| `LD_FRONTIER_AXES` | `tail,hetero,nonlinear` | Axes to sweep |
| `LD_OUTER_REPS` | `15` | Audit replications; `0` disables audit mode |
| `LD_OUTER_CELLS` | `signed:wellspec,signed:student_t` | Audit cells as `regret:dgp` pairs |
| `LD_OUTER_KEEP` | `0` | `1` keeps per-replication scratch folders |
| `LD_REGRET_MODE` | unset | `positive` or `signed`; grid mode runs both when unset |
| `LD_DGP_MODE` | `all` | Grid-mode data process: `wellspec`, `student_t`, `heteroskedastic`, `nonlinear`, or `all` |
| `LD_DGP_TAIL_DF`, `LD_DGP_HET_THETA`, `LD_DGP_NL_COEF` | `Inf`, `0`, `0` | Manual severity knobs |
| `LD_N_CAL`, `LD_N_TUNE`, `LD_N_TEST` | per mode | Inner paths per cell (grid 20/8/100, audit 12/8/20, frontier 8/6/10) |
| `LD_OUTPUT_DIR` | home | Output root |
| `LD_RESUME` | `1` | `0` disables checkpoint resumption |

## Mapping the manuscript to outputs

| Manuscript | Source |
| --- | --- |
| Table 1, Section 5.1 | `table1_primary_q75_v13d.csv` in the signed well-specified grid folder |
| Section 5.2 sensitivities | `relative_summary_v13d.csv` in the signed and positive grid folders |
| Table 2, Section 5.3 | `policy_summary_v13d.csv` in the positive grid folder (also committed under `outputs/`) |
| Table 3, Section 5.4 | `outer_mc_summary_v13d.csv` in the two audit folders |
| Table 4, bands, crossings, diagnostic (Section 6) | `frontier_summary_v13d.csv` and `frontier_crossings_v13d.csv` in the frontier folder |
| Appendix A figures | Written by grid mode when `ggplot2` is installed |

## Verification

```
Rscript verify_paper_numbers.R
```

Recomputes every number reported in the manuscript, 183 checks in total, and
prints PASS or FAIL per value. For each dataset it looks first under
`replication/` in this repository, then under `LD_OUTPUT_ROOT` (default: the
home directory), and for the positive primary run it also falls back to the
committed `outputs/` directory. Missing folders skip their group with an
explanation rather than failing, and the script exits with status 1 if any
available check fails, so it can gate a release.

## Determinism and resumption

Seed bases are fixed per stage, and each audit or frontier replication uses an
offset seed stream, so every study is reproducible bit for bit. The
well-specified baseline preserves the random-number stream of the original
committed run exactly. All long stages checkpoint to their output folder and
resume where they stopped.
