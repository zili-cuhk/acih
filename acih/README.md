# acih — code for "Active-Coordinate Inference for Sparse Additive Hazards Regression under Informative Subsampling"

R package (base R only: `stats`, `utils`, `graphics`, `grDevices`) implementing the method, the simulation
study of Section 4 (Studies 1–4), the data-analysis protocol of Section 5 and the plasmode study. Written for R 4.5.1 on an LSF cluster.

## Package layout

| Path | Content |
| --- | --- |
| `R/00_utils.R` | helpers (AR(1)/block covariates, ridge solve, Wilson interval) |
| `R/01_ah_core.R` | Lin–Ying estimating equations with risk-set centring: blockwise ambient `b`, `V` for screening (`ah_bV`); per-subject quantities on a selected set (`ah_setup`), `q_i(beta)`, Breslow residual `xi_i(beta)`; split-full refit with martingale sandwich |
| `R/02_ipw.R` | with-replacement allocation rules (uniform, score norm, Zuo et al., inference-oriented with contrast matrix C), event stratification and stabilisation, inverse-CDF draws, IPW refit and total sandwich |
| `R/03_datagen.R` | additive-hazards data generation (bounded covariates, nonnegative hazard, constant or linear baseline, time-varying effect), censoring (uniform, covariate-dependent exponential, Weibull, shared frailty, none), calibration, least-false targets |
| `R/04_screen.R` | GSELO continuation path, joint EBIC (BIC term based on the number of screening-split events, size cap 20 in the simulations and 40 in the data analyses), marginal pre-screen, fixed-budget ranking |
| `R/05_simrep.R` | one replicate of Study 1 (`run_rep_known`), Studies 2–3 (`run_rep_pipeline`) and the marginal route (`run_rep_marginal`) |
| `R/06_realdata.R` | data-analysis protocol: one split per task (`run_app_task`), summaries and figures (`summarize_apps`) |
| `R/07_extensions.R` | Poisson design, including the event-complete design, HT variance; the case-cohort comparator of Kulich and Lin (2000) |
| `R/08_designs.R` | the designs of Studies 1–3 (`design_table`, `prepare_designs`) |
| `R/09_run.R` | task table, task runner, timing, run directory, smoke test |
| `R/10_summarize.R` | tables (CSV, LaTeX) and figures |
| `R/12_figures.R` | figure style and drawing; `make_figures_from_csv()` (`Rscript acih.R figures`) redraws all figures from the summary CSV files |
| `R/11_plasmode.R` | plasmode study on the processed cohorts (support-known refits with real covariates) |
| `R/13_tables.R` | `make_tables_from_csv()` (`Rscript acih.R tables`) writes every table of the manuscript and supplement from the summary CSV files |
| `inst/scripts/acih.R` | command-line entry point |
| `inst/hpc/` | LSF scripts (queue names and limits may need local adjustment) |
| `inst/app/` | MSK-MET (clinical and genomic) and NHIS-LMF preprocessing |
| `tests/run_tests.R` | exactness of the blockwise ambient quantities, score identities, conditional unbiasedness and exact variance of the with-replacement and Poisson designs, smoke test |

## Simulation study (Section 4)

* **Study 1** (refit stage, selected set `S0 + 4` noise covariates, 1,000 replicates): `m` in {10k, 25k, 50k} x
  {50% uniform censoring, 50% covariate-dependent censoring with rho = 0.75, rare events with 90% censoring};
  4% events without censoring (`m` = 25k, 50k); linear baseline; ten active coefficients; binary covariates.
  Methods: SplitFull, UIPW, SIPW, ZIPW, IIPW, CIPW, PUIPW, PIIPW, CaseCohort, the Omega-only ablation of IIPW.
* **Study 2** (complete pipeline, 500 replicates): the seven primary designs of the original draft
  (`n` in {50k, 100k, 200k} x `p` in {500, 1000}, and `n` = 100k, `p` = 500, rho = 0.75 with covariate-dependent
  censoring); SplitFull, UIPW, SIPW, IIPW, PUIPW, PIIPW, CaseCohort at `r` in {1000, 5000}, OracleFull, and Naive.
* **Study 3** (500 replicates): screening challenges (block correlation, bounded t5 with Weibull censoring, weak
  signals, ten active coefficients, rare events, 4% events, binary covariates); the `p >= n_scr` marginal route;
  misspecification (time-varying coefficient, dependent censoring through a shared frailty) with least-false
  targets.

* **Study 4** (complete pipeline, 200 replicates): biobank scale, `n` = 500,000, `p` = 500, about 4% and 2% events
  (no censoring before tau), `r` in {10000, 20000, 50000}.
* **Plasmode study** (500 replicates per cohort): real covariates of the three processed cohorts, true support =
  variables selected in at least 95 of 100 splits of the data analysis, refitted set = primary-split selected set;
  run with `bash hpc/submit_plasmode.sh` after the data analysis.
* Every IPW/Poisson fit also records the subsampling-only (conditional) variance of the existing subsampling
  literature (`ZIPW_GammaOnly` = the Zuo-type allocation with the conditional variance of Zuo et al. 2021, `IIPW_GammaOnly`, `PIIPW_GammaOnly`) and the
  coverage of the full-inference-sample estimate (`CP_full`).

## Running on an LSF cluster

```bash
## 1. install (once)
module load R/4.5.1
mkdir -p ~/Rlib
R CMD INSTALL -l ~/Rlib acih
export R_LIBS=~/Rlib${R_LIBS:+:$R_LIBS}

## 2. run directory
Rscript -e 'acih::setup_run_dir("acih_run")'
cd acih_run
Rscript acih.R test                     # smoke test

## 3. calibrate the designs, then time the heaviest designs (recommended before the full array)
bsub < hpc/00_prepare.lsf               # writes out/designs.rds
bsub -w "done(acih_prepare)" < hpc/04_time.lsf
## look at logs/time.*.err: seconds per replicate and hours per task

## 4. the simulation array and the summary
bash hpc/submit_all.sh                  # needs R on the PATH: run 'module load R/4.5.1' first
```

`Rscript acih.R list` prints the task table (2,350 tasks; the `%400` in the array limits concurrent jobs).
Results: `out/<study>/<design>/reps_*.rds`; summaries in `results/` (`coordinate_level.csv`,
`coefficient_averaged.csv`, `screening.csv`, `marginal_route.csv`, `fig_S1_calibration.pdf`,
`fig_S1_efficiency.pdf`, `fig_failure_mechanisms.pdf`, `fig_S2_calibration.pdf`, `fig_calibration_summary.pdf`, `table_S1.tex`, `table_S2.tex`, `table_screening.tex`).
Rerunning a task skips replicates whose output file exists, so failed tasks can simply be resubmitted.

## Data applications (Section 5)

Raw files go under `app/data/` in the run directory:

| Path | Source |
| --- | --- |
| `app/data/msk_met_2021/data_clinical_patient.txt`, `data_clinical_sample.txt` | cBioPortal datahub, `public/msk_met_2021` (Git LFS) |
| `app/data/msk_met_2021/data_mutations.txt` | mutation calls of `msk_met_2021` from the cBioPortal web API: `Rscript app/msk_met_fetch_mutations.R` (needs internet, `curl` and the `jsonlite` package) writes the three columns used (`Hugo_Symbol`, `Tumor_Sample_Barcode`, `Variant_Classification`; 227,931 calls). The DataHub file (Git LFS) can be used instead |
| `app/data/gene_panels/impact341.txt`, `impact410.txt`, `impact468.txt` | datahub `reference_data/gene_panels/data_gene_panel_impact*.txt` |
| `app/data/nhis/samadult_<y>.csv`, `person_<y>.csv`, `NHIS_<y>_MORT_2019_PUBLIC.dat`, y = 2016–2018 | NCHS NHIS public-use files and Linked Mortality Files |

```bash
bash hpc/submit_apps.sh      # preprocessing (waits), array over cohorts x 100 splits, summary
```

The preprocessing log must show 24,561 patients, 10,095 deaths and 10,025 deaths within 60 months for MSK-MET,
and 83,676 adults, 2,739 deaths and 2,545 deaths within three years for NHIS-LMF. Cohorts: `msk_met` (20
clinical candidates), `msk_met_genomic` (the clinical candidates plus mutation indicators of genes on all
IMPACT panels mutated in at least 1% of patients) and `nhis_lmf`. Split 1 is the prespecified primary split
(seed 20260701) and also carries the fixed-budget (K = 5, 10, 15) and survey-weighted sensitivity analyses.
Results: `app/out/<cohort>/split_*.rds`; tables and figures in `app/results/`.
`app_summarize` also reruns the screen of the genomic analysis with the simulation cap of 20 (`app_cap_sensitivity()`,
also available as `Rscript acih.R app_cap20`), which gives `msk_met_genomic_cap20_*.csv`.

The GSELO screen is computed on the time scale on which the crude event rate of the screening data on
`[0, tau]` equals one (`time_scale = "rate"` in `screen_gselo`), so the selected set does not depend on the unit
of time; refits and inference are equivariant under rescaling of time.

## Notes

* `ZIPW` implements a Zuo-type allocation, proportional to an upper bound of the score norm as in Zuo et al. (2021);
  it is not their complete procedure.
* `CaseCohort` is the estimator and variance estimator of Kulich and Lin (2000, Biometrika 87, 73-87, eqs. 3, 6, 7)
  for a Bernoulli case-cohort design. It is computed on exactly the same event-complete subsample as `PUIPW`
  (all cases, controls with probability `(r - n1)/(m - n1)`), but, as in a case-cohort study, its risk-set
  averages use only the sampled records; it is not defined (reported as missing) when the number of events
  exceeds `r`.
* `RIPW` (residual-score IPW) is exploratory and not part of the manuscript.
* A fixed `1e-8` ridge is used in every linear solve, matching the manuscript's numerical safeguard.

## Method names

Package identifiers and the compact names used in the paper, figures, and tables:

| Identifier | Name | Description |
|---|---|---|
| `SplitFull` | FULL | Full inference-split refit (no subsampling) |
| `OracleFull` | ORACLE | Full inference-split refit on the true support |
| `Naive` | NAIVE | Screening and refitting on the same full sample |
| `UIPW` | UNIF | Uniform sampling with replacement |
| `SIPW` | SNORM | Score-norm sampling with replacement |
| `ZIPW` | ZUO | Zuo-type allocation (Zuo et al., 2021) |
| `IIPW` | OPT | A-optimal sampling with replacement |
| `CIPW` | OPT-C | A-optimal sampling targeted at a contrast `C_A` (first active coordinate) |
| `PUIPW` | P-UNIF | Uniform Poisson sampling |
| `PIIPW` | P-OPT | A-optimal Poisson sampling (recommended) |
| `CaseCohort` | KL-CC | Kulich-Lin case-cohort on the P-UNIF subsample |

Suffixes `_GammaOnly` / `_OmegaOnly` mark the Gamma-only (conditional) and Omega-only variance versions; no suffix is the total variance.
