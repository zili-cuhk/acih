## =====================================================================
## 08_designs.R -- the simulation study of the manuscript (Section 4).
##  S1  Support-known refit stage (A = S0 + 4 noise covariates): calibration
##      and efficiency of every allocation rule across (m, r, censoring,
##      baseline, dimension of S0, covariate type); variance ablation;
##      contrast-targeted allocation.  1000 replications.
##  S2  Full pipeline with GSELO screening (the seven primary designs of the
##      manuscript), conditional and pipeline coverage, OracleFull, Naive
##      (no splitting) reference.  500 replications.
##  S3  Stress designs: screening challenges (block, bounded t5, weak
##      signals, s = 10, rare events, binary covariates), the p >= n_scr
##      marginal route, and model misspecification (time-varying effect,
##      dependent censoring) with least-false targets.  500 replications.
## =====================================================================

design_table <- function() {
  L <- list()
  add <- function(id, study, des_args, run, R, extra = list())
    L[[id]] <<- list(id = id, study = study, des_args = des_args, run = run, R = R, extra = extra)
  ## ---- S1 ----
  for (m in c(10000, 25000, 50000)) {
    rv <- unique(c(1000, 5000, m / 2))
    add(sprintf("S1_m%dk_U50", m / 1000), "S1", list(p = 8, rho = 0.5, cens = "uniform", cens_rate = 0.5), "known", 1000, list(m = m, r = rv))
    add(sprintf("S1_m%dk_CD50", m / 1000), "S1", list(p = 8, rho = 0.75, cens = "covexp", cens_rate = 0.5), "known", 1000, list(m = m, r = rv))
    ## rare events (about 10% observed events): short follow-up tau = 0.03,
    ## uniform censoring calibrated to 90%
    add(sprintf("S1_m%dk_U90", m / 1000), "S1", list(p = 8, rho = 0.5, cens = "uniform", cens_rate = 0.9, tau = 0.03), "known", 1000, list(m = m, r = rv))
  }
  rv25 <- c(1000, 5000, 12500)
  ## the linear baseline and the ten-coefficient design have larger hazards; their
  ## horizons are shortened so that about 10-15% of subjects remain at risk at tau
  ## (tau = 0.3 would leave fewer than the 5% required by the risk-set gate)
  add("S1_m25k_LIN", "S1", list(p = 8, rho = 0.5, cens = "uniform", cens_rate = 0.5, base_type = "linear", slope_frac = 3, tau = 0.2), "known", 1000, list(m = 25000, r = rv25))
  add("S1_m25k_S10", "S1", list(p = 14, s_design = "s10", rho = 0.5, cens = "uniform", cens_rate = 0.5, tau = 0.12), "known", 1000, list(m = 25000, r = rv25))
  ## event-complete Poisson design is most relevant with rare events: about
  ## 4% events (no censoring before tau = 0.01, administrative stop at tau)
  for (m in c(25000, 50000))
    add(sprintf("S1_m%dk_RARE4", m / 1000), "S1", list(p = 8, rho = 0.5, cens = "none", tau = 0.01), "known", 1000,
        list(m = m, r = if (m == 25000) c(1500, 2500, 5000) else c(3000, 5000, 10000)))
  add("S1_m25k_BIN", "S1", list(p = 8, rho = 0.5, cens = "uniform", cens_rate = 0.5, binary_every = 2), "known", 1000, list(m = 25000, r = rv25))
  ## ---- S2 ----
  for (n in c(50000, 100000, 200000)) for (p in c(500, 1000))
    add(sprintf("S2_n%dk_p%d_U50", n / 1000, p), "S2", list(n = n, p = p, rho = 0.5, cens = "uniform", cens_rate = 0.5),
        "pipeline", 500, list(naive = (n <= 100000 && p == 500)))
  add("S2_n100k_p500_CD50", "S2", list(n = 100000, p = 500, rho = 0.75, cens = "covexp", cens_rate = 0.5), "pipeline", 500, list(naive = TRUE))
  ## ---- S3: screening stress (pipeline) ----
  b <- list(n = 100000, p = 500)
  add("S3_block", "S3", c(b, list(xdesign = "block", rho = 0.75, cens = "uniform", cens_rate = 0.5)), "pipeline", 500)
  add("S3_t5", "S3", c(b, list(xdesign = "t5", rho = 0.75, cens = "weibull", cens_rate = 0.5)), "pipeline", 500)
  add("S3_weak", "S3", c(b, list(s_design = "weak4", rho = 0.5, cens = "uniform", cens_rate = 0.5)), "pipeline", 500)
  add("S3_s10", "S3", c(b, list(s_design = "s10", rho = 0.5, cens = "uniform", cens_rate = 0.5, tau = 0.12)), "pipeline", 500)
  add("S3_U90", "S3", c(b, list(rho = 0.5, cens = "uniform", cens_rate = 0.9, tau = 0.03)), "pipeline", 500)
  add("S3_RARE4", "S3", c(b, list(rho = 0.5, cens = "none", tau = 0.01)), "pipeline", 500, list(r = c(5000, 10000)))
  add("S3_bin", "S3", c(b, list(rho = 0.5, cens = "uniform", cens_rate = 0.5, binary_every = 5)), "pipeline", 500)
  ## ---- S3: p >= n_scr marginal route ----
  pn <- list(n = 1000, p = 2000, cens = "uniform", cens_rate = 0.25)
  add("S3_pn_indep_alt", "S3", c(pn, list(s_design = "alt4", rho = 0)), "marginal", 500)
  add("S3_pn_ar50_reinf", "S3", c(pn, list(s_design = "reinf4", rho = 0.5)), "marginal", 500)
  add("S3_pn_ar85_alt", "S3", c(pn, list(s_design = "alt4", rho = 0.85)), "marginal", 500)
  add("S3_pn_block_reinf", "S3", c(pn, list(s_design = "reinf4", xdesign = "block", rho = 0.75)), "marginal", 500)
  add("S3_pn_t5_reinf", "S3", c(pn, list(s_design = "reinf4", xdesign = "t5", rho = 0.5)), "marginal", 500)
  ## ---- S4: biobank-scale rare events (pipeline) ----
  ## n = 500,000 and p = 500; no censoring before tau, about 4% (tau = 0.01)
  ## and 2% (tau = 0.005) events, so that the screening split carries about
  ## 9,300 and 4,900 events; the event-complete design applies at every r
  for (k in list(c(4, 0.01), c(2, 0.005)))
    add(sprintf("S4_bb%d", k[1]), "S4", list(n = 500000, p = 500, rho = 0.5, cens = "none", tau = k[2]),
        "pipeline", 200, list(r = c(10000, 20000, 50000)))
  ## ---- S3: misspecification (support known, least-false target) ----
  add("S3_tv", "S3", list(p = 8, rho = 0.5, cens = "uniform", cens_rate = 0.5,
                          tv = list(j = 1, slope = 4, center = 0.15)), "known", 500,
      list(m = 25000, r = c(1000, 5000), lf_target = TRUE))
  add("S3_frailty", "S3", list(p = 8, rho = 0.5, cens = "frailty", cens_rate = 0.5), "known", 500,
      list(m = 25000, r = c(1000, 5000), lf_target = TRUE))
  L
}

## Build, calibrate and (where needed) compute least-false targets once
prepare_designs <- function(dir = ".", n_cal = 2e5, n_big = 1e6, verbose = TRUE) {
  dir.create(file.path(dir, "out"), showWarnings = FALSE, recursive = TRUE)
  file <- file.path(dir, "out", "designs.rds")
  tab <- design_table()
  out <- if (file.exists(file)) readRDS(file) else list()
  for (x in tab) {
    if (!is.null(out[[x$id]])) next
    args <- x$des_args; if (is.null(args$n)) args$n <- x$extra$m %||% 1e5
    des <- do.call(make_design, c(args, list(label = x$id)))
    des <- calibrate_design(des, n_cal = n_cal)
    if (isTRUE(x$extra$lf_target)) {
      A <- c(des$S0, setdiff(seq_len(des$p), des$S0)[1:4])
      des$lf_target <- population_target(des, A, n_big = n_big)[seq_along(des$S0)]
    }
    des$study <- x$study; des$run <- x$run; des$R <- x$R; des$extra <- x$extra
    out[[x$id]] <- des
    if (verbose) message(sprintf("%-22s cens=%.3f risk(tau)=%.3f events(tau)=%.3f", x$id,
                                 des$achieved_cens, des$risk_frac_tau, des$event_frac_tau))
    saveRDS(out, file)
  }
  invisible(out)
}
