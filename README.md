# acih: active-coordinate inference for sparse additive hazards regression under informative subsampling

Code and result summaries for the manuscript

> Zili Liu, Juan Liu, and Hong Wang. *Active-Coordinate Inference for Sparse Additive Hazards Regression under
> Informative Subsampling.* Submitted to *Statistics in Medicine* (2026).

The method screens covariates on one half of the data and then refits the additive hazards model of Lin and
Ying (1994) on a weighted subsample of the other half. The confidence intervals for the true coefficients use a
total sandwich covariance that adds the subsampling variation to the variation of the inference split. Designs
include sampling with replacement, Poisson sampling, and an event-complete Poisson design that keeps every
observed event.

## Contents

| Path | Content |
| --- | --- |
| `acih/` | R package (base R only) with the method, the simulation studies, the data-analysis protocol, the plasmode study, LSF cluster scripts and tests; see [`acih/README.md`](acih/README.md) |
| `results/simulation/` | summaries of Studies 1–4: coordinate-level and coefficient-averaged results, screening diagnostics, the marginal pre-screen and the comparison of the two screening criteria |
| `results/applications/` | summaries of the MSK-MET and NHIS-LMF analyses (selected sets, estimates, efficiency, sensitivity analyses, Kaplan–Meier summaries) and of the plasmode study |
| `figures/` | the figures of the manuscript and its supplement, drawn from the files in `results/` |
| `tables/` | the tables of the manuscript and its supplement (LaTeX `tabular` environments), computed from the files in `results/` |

No individual-level data are included.

## Installation

```bash
git clone https://github.com/zili-cuhk/acih.git
R CMD INSTALL acih/acih
```

The package needs R (>= 4.3.0) and no other packages. `Rscript acih/tests/run_tests.R` checks the
estimating equations, the conditional unbiasedness and exact variances of the subsampling designs, and runs one
small replicate of each study type (about a minute).

## Quick example

```r
library(acih)
des <- calibrate_design(make_design(n = 20000, p = 50, cens = "uniform", cens_rate = 0.5), n_cal = 2e4)
d <- gen_data(des, seed = 1)
## 50/50 split: GSELO screen on one half, refits on the other half at budgets r = 1000 and 4000
sp <- primary_split(nrow(d$X), 1)
fit <- analyze_split(d$time, d$status, d$X, des$tau, sp, "gselo", r_vec = c(1000, 4000))
fit$table   # estimates and total-sandwich standard errors for every design and budget
```

## Reproducing the figures and tables from the summaries

```r
library(acih)
make_figures_from_csv(sim_dir = "results/simulation", app_dir = "results/applications", out = "figures")
make_tables_from_csv(sim_dir = "results/simulation", app_dir = "results/applications", out = "tables")
```

This redraws every figure and rewrites every table in a few seconds without refitting.

| Figure | File |
| --- | --- |
| Figure 1 | `fig_calibration_summary.pdf` |
| Figure 2 | `fig_failure_mechanisms.pdf` |
| Figure 3 | `fig_S1_efficiency.pdf` |
| Figure 4 | `fig_real_data_efficiency.pdf` |
| Figure 5 | `fig_plasmode.pdf` |
| Web Figures 1–3 | `fig_S1_calibration.pdf`, `fig_S2_calibration.pdf`, `fig_real_data_km.pdf` |

| Table | File |
| --- | --- |
| Table 1 | `tab_sim_design.tex` |
| Table 2 | `tab_real_data_primary.tex` |
| Table 3 | `tab_real_data_efficiency.tex` |
| Web Tables 1–6 (Study 1) | `tab_supp_s1_cp.tex`, `tab_supp_s1_ratio.tex`, `tab_supp_s1_ail.tex`, `tab_supp_variance_ablation.tex`, `tab_supp_contrast.tex`, `tab_supp_gamma_only.tex` |
| Web Tables 7–8 (Study 2) | `tab_supp_s2.tex` |
| Web Tables 9–13 (Study 3) | `tab_supp_screening_challenges.tex`, `tab_supp_ebic_compare.tex`, `tab_supp_misspec.tex`, `tab_supp_marginal_prescreen_pgeqn.tex` |
| Web Table 14 (Study 4) | `tab_supp_s4.tex` |
| Web Tables 15–22 (data analyses) | `tab_supp_real_screening.tex`, `tab_supp_primary_msk_met.tex`, `tab_supp_primary_nhis_lmf.tex`, `tab_supp_cap_genomic.tex`, `tab_supp_budget_*.tex`, `tab_supp_nhis_survey_sensitivity.tex` |
| Web Tables 23–25 (plasmode) | `tab_supp_plasmode.tex` |

Web Tables 7–8, 9–10 and 23–25 are single tables split over pages in the supplement.

## Reproducing the full study

The simulation study (about 2,350 cluster tasks), the data analysis (100 splits per cohort) and the plasmode
study (500 replications per cohort) were run with R 4.5.1 on an LSF cluster. [`acih/README.md`](acih/README.md)
describes the run directory, the job scripts in `acih/inst/hpc/`, the seeds, and the expected output. In short:

```bash
Rscript -e 'acih::setup_run_dir("acih_run")'
cd acih_run
bash hpc/submit_all.sh        # Studies 1-4 and their summaries
bash hpc/submit_apps.sh       # data analyses (raw files under app/data/, see below)
bash hpc/submit_plasmode.sh   # plasmode study, after the data analyses
bash hpc/submit_screencmp.sh  # comparison of the two screening criteria
Rscript acih.R figures        # figures from the summaries
Rscript acih.R tables         # tables from the summaries
```

## Data

* **MSK-MET** (`msk_met_2021`): clinical files from the [cBioPortal DataHub](https://github.com/cBioPortal/datahub)
  and mutation calls from the DataHub or the cBioPortal web API; preprocessing in `acih/inst/app/msk_met_prepare.R`.
* **NHIS-LMF**: 2016–2018 National Health Interview Survey sample-adult and person files and the 2019 public-use
  Linked Mortality Files from the National Center for Health Statistics; preprocessing in
  `acih/inst/app/nhis_lmf_prepare.R`.

Both data sources are public and de-identified. Users must obtain them from the providers under their terms.

## License

MIT; see [`LICENSE`](LICENSE).
