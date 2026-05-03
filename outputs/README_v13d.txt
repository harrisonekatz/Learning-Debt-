Learning Debt simulation outputs, v13d

Run mode: paper
T_PERIODS: 200
Calibration paths per scenario cell: 20
Tuning paths per scenario cell: 8
Test paths per scenario cell: 100

Design summary:
- Warm-started deployed and shadow posteriors.
- Separate update, monitor, and evaluation batches.
- Period-t actions affect period t+1 and later.
- Objective is lambda times retrain count plus accumulated positive evaluation score gap.
- Exact KL between NIG shadow and deployed posteriors is the debt signal.
- Debt is age-adjusted using stable no-shift calibration paths.
- Utility-calibrated debt and hybrid policies are tuned on calibration paths and evaluated on held-out paths.
- Calendar, CUSUM, alarm, debt-threshold, always-retrain, and never-retrain baselines are included.
- v13d expands the grids because v13c placed some baselines at boundary values.
- Score-unit sensitivity is run for median, q75, and mean positive predictive-regret units.
- The script checkpoints calibration, tuning, and held-out test paths and can resume after interruption.

Expanded grids:
- Calendar intervals: 1, 2, 3, 5, 10, 20, 40, 80
- CUSUM thresholds: 0.01, 0.025, 0.05, 0.1, 0.2, 0.5, 1, 2, 4, 8, 16
- Alarm quantiles: 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975
- Utility scales: 0.25, 0.5, 0.75, 1, 1.5, 2, 4, 8
- Debt-threshold quantiles: 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99

Main files:
- all_runs_v13d.csv
- policy_summary_v13d.csv
- relative_summary_v13d.csv
- cell_win_summary_v13d.csv
- replicate_win_summary_v13d.csv
- score_unit_sensitivity_summary_v13d.csv
- tuning_summary_v13d.csv
- tuned_parameters_v13d.csv
- boundary_summary_v13d.csv
- calibration_rows_v13d.csv
- score_units_v13d.csv
- session_info_v13d.txt
