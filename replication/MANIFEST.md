# Replication outputs manifest

Copy each output folder from your output root (by default the home directory
on the machine that ran the studies) into this `replication/` directory,
keeping the folder name unchanged. Only the files listed are required by
`verify_paper_numbers.R`; the per-replication files are included in the
submission package for completeness.

| Folder | Required files | Status |
| --- | --- | --- |
| `learning_debt_sim_outputs_v13d_paper_signed/` | `score_units_v13d.csv`, `table1_primary_q75_v13d.csv`, `relative_summary_v13d.csv` | copy from run machine |
| `learning_debt_sim_outputs_v13d_paper_positive/` | `relative_summary_v13d.csv`, `policy_summary_v13d.csv` | already available via the committed `outputs/` directory; copy here only if you want everything in one place |
| `learning_debt_sim_outputs_v13d_paper_signed_outermc/` | `outer_mc_summary_v13d.csv`, `outer_mc_per_rep_v13d.csv` | copy from run machine |
| `learning_debt_sim_outputs_v13d_paper_signed_student_t_outermc/` | `outer_mc_summary_v13d.csv`, `outer_mc_per_rep_v13d.csv` | INCLUDED in this package |
| `learning_debt_sim_outputs_v13d_paper_signed_frontier/` | `frontier_summary_v13d.csv`, `frontier_crossings_v13d.csv`, `frontier_per_rep_v13d.csv` | copy from run machine |

Optionally also copy the remaining grid folders (the positive and signed
student_t, heteroskedastic, and nonlinear cells) if you want the
misspecification variants quoted in review correspondence to travel with the
repository.

After copying, run `Rscript verify_paper_numbers.R` from the repository root.
All 183 checks should pass with no skipped groups.
