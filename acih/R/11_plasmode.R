## =====================================================================
## File 11_plasmode.R -- plasmode simulation built on a real cohort.
## The covariate matrix of the cohort is kept fixed; outcomes are generated
## from an additive hazards model whose coefficients are the split-free
## fit of a reference support (the primary-split selected set of the data
## analysis) on the whole cohort, with all other coefficients zero:
##   lambda(t | X) = lambda0 + X^T beta,
## lambda0 = max(crude event rate, -min_i X_i^T beta + floor) so that every
## record has a positive hazard; uniform censoring calibrated to the observed
## proportion of events within tau; administrative censoring at tau.
## Each replicate draws new outcomes, a new 50/50 split, screens with GSELO,
## and refits the selected set (split-full, with-replacement and Poisson
## designs, case-cohort, conditional-variance comparators, OracleFull).
## =====================================================================

plasmode_design <- function(dat, covs, support, tau, floor_rate = NULL, n_cal = 2e5, seed = 1) {
  X <- as.matrix(dat[, covs, drop = FALSE]); colnames(X) <- covs
  S0 <- match(support, covs)
  if (anyNA(S0)) stop("support not in the candidate set")
  st <- ah_setup(dat$time, dat$status, X[, S0, drop = FALSE], tau)
  beta <- numeric(length(covs)); beta[S0] <- ah_fit_full(st)$beta
  lp <- drop(X %*% beta)
  ev_obs <- mean(dat$status == 1 & dat$time <= tau)
  crude <- sum(dat$status == 1 & dat$time <= tau) / sum(pmin(dat$time, tau))
  floor_rate <- floor_rate %||% (0.05 * crude)
  lambda0 <- max(crude, -min(lp) + floor_rate)
  ## calibrate uniform censoring C ~ U(0, c) to the observed event proportion within tau
  set.seed(seed); i <- sample.int(nrow(X), n_cal, replace = TRUE)
  h <- lambda0 + lp[i]; Tt <- stats::rexp(n_cal, h); V <- stats::runif(n_cal)
  evp <- function(c) mean(Tt <= pmin(V * c, tau))
  c_cens <- if (evp(1e6 * tau) <= ev_obs) Inf else
    exp(stats::uniroot(function(lc) evp(exp(lc)) - ev_obs, c(log(tau * 1e-3), log(tau * 1e6)))$root)
  structure(list(X = X, beta = beta, S0 = S0, lambda0 = lambda0, c_cens = c_cens, tau = tau,
                 lp_range = range(lp), crude = crude, ev_obs = ev_obs, ev_sim = evp(c_cens)),
            class = "plasmode_design")
}

plasmode_data <- function(pd) {
  n <- nrow(pd$X); h <- pd$lambda0 + drop(pd$X %*% pd$beta)
  Tt <- stats::rexp(n, h)
  C <- if (is.finite(pd$c_cens)) stats::runif(n, 0, pd$c_cens) else rep(Inf, n)
  time <- pmin(Tt, C, pd$tau); status <- as.integer(Tt <= pmin(C, pd$tau))
  list(time = time, status = status)
}

run_rep_plasmode <- function(pd, seed, r_vec, cap = 40,
                             methods = ipw_methods_default()[c("UIPW", "ZIPW", "IIPW", "PUIPW", "PIIPW")]) {
  set.seed(seed)
  d <- plasmode_data(pd); X <- pd$X; n <- nrow(X); S0 <- pd$S0
  scr <- sort(sample.int(n, floor(n / 2))); inf <- setdiff(seq_len(n), scr)
  sc <- screen_gselo(d$time[scr], d$status[scr], X[scr, , drop = FALSE], pd$tau, cap = cap)
  A <- sc$A
  sdiag <- c(size = length(A), sure = all(S0 %in% A), exact = setequal(A, S0), n_false = length(setdiff(A, S0)),
             t_screen = sc$time, events_inf = sum(d$status[inf]))
  out <- list()
  if (length(A)) {
    pil <- pilot_fit(d$time[scr], d$status[scr], X[scr, A, drop = FALSE], pd$tau)
    st <- ah_setup(d$time[inf], d$status[inf], X[inf, A, drop = FALSE], pd$tau)
    f <- ah_fit_full(st)
    out[[1]] <- .rec("SplitFull", NA, f$beta, f$se, A, S0)
    U <- stats::runif(max(r_vec)); Up <- stats::runif(st$n)
    for (mn in names(methods)) {
      mm <- methods[[mn]]
      if (identical(mm$design, "poisson")) { out <- c(out, .run_poisson(mn, mm, st, pil, r_vec, Up, A, S0, 0)); next }
      pi <- ipw_probs(st, pil, mm$rule, score = mm$score)
      for (r in r_vec) {
        g <- ah_fit_ipw(st, pi, draw_ipw(pi, U[seq_len(r)]), score = mm$score)
        out[[length(out) + 1]] <- .rec(mn, r, g$beta, g$se, A, S0, list(gamma_frac = g$gamma_frac))
        if (mn == "IIPW") out[[length(out) + 1]] <- .rec("IIPW_OmegaOnly", r, g$beta, sqrt(diag(g$Sigma_omega_only)), A, S0)
        if (mn %in% c("IIPW", "ZIPW"))
          out[[length(out) + 1]] <- .rec(paste0(mn, "_GammaOnly"), r, g$beta, sqrt(diag(g$Sigma_gamma_only)), A, S0)
      }
    }
    fb <- rep(NA_real_, length(S0)); pos <- match(S0, A); fb[!is.na(pos)] <- f$beta[pos[!is.na(pos)]]
  } else {
    fb <- rep(NA_real_, length(S0))
    out[[1]] <- .rec("SplitFull", NA, numeric(0), numeric(0), integer(0), S0)
    for (mn in names(methods)) for (r in r_vec) out[[length(out) + 1]] <- .rec(mn, r, numeric(0), numeric(0), integer(0), S0)
  }
  st0 <- ah_setup(d$time[inf], d$status[inf], X[inf, S0, drop = FALSE], pd$tau); f0 <- ah_fit_full(st0)
  out[[length(out) + 1]] <- .rec("OracleFull", NA, f0$beta, f0$se, S0, S0)
  res <- do.call(rbind, out); res$true <- pd$beta[S0][res$coord]; res$full <- fb[res$coord]
  list(gate_fail = FALSE, res = res, screen = sdiag)
}

## Support-known plasmode replicate: new outcomes, new 50/50 split, pilot on
## one half, refit of the fixed reference set `A` (names) on the other half.
run_rep_plasmode_known <- function(pd, A, seed, r_vec,
                                   methods = ipw_methods_default()[c("UIPW", "SIPW", "ZIPW", "IIPW", "PUIPW", "PIIPW")]) {
  set.seed(seed)
  d <- plasmode_data(pd); X <- pd$X; n <- nrow(X); S0 <- pd$S0
  Ai <- match(A, colnames(X)); if (anyNA(Ai) || !all(S0 %in% Ai)) stop("A must contain the true support")
  scr <- sort(sample.int(n, floor(n / 2))); inf <- setdiff(seq_len(n), scr)
  pil <- pilot_fit(d$time[scr], d$status[scr], X[scr, Ai, drop = FALSE], pd$tau)
  st <- ah_setup(d$time[inf], d$status[inf], X[inf, Ai, drop = FALSE], pd$tau)
  f <- ah_fit_full(st)
  out <- list(.rec("SplitFull", NA, f$beta, f$se, Ai, S0))
  U <- stats::runif(max(r_vec)); Up <- stats::runif(st$n)
  for (mn in names(methods)) {
    mm <- methods[[mn]]
    if (identical(mm$design, "poisson")) { out <- c(out, .run_poisson(mn, mm, st, pil, r_vec, Up, Ai, S0, 0)); next }
    pi <- ipw_probs(st, pil, mm$rule, score = mm$score)
    for (r in r_vec) {
      g <- ah_fit_ipw(st, pi, draw_ipw(pi, U[seq_len(r)]), score = mm$score)
      out[[length(out) + 1]] <- .rec(mn, r, g$beta, g$se, Ai, S0, list(gamma_frac = g$gamma_frac))
      if (mn == "IIPW") out[[length(out) + 1]] <- .rec("IIPW_OmegaOnly", r, g$beta, sqrt(diag(g$Sigma_omega_only)), Ai, S0)
      if (mn %in% c("IIPW", "ZIPW"))
        out[[length(out) + 1]] <- .rec(paste0(mn, "_GammaOnly"), r, g$beta, sqrt(diag(g$Sigma_gamma_only)), Ai, S0)
    }
  }
  res <- do.call(rbind, out); res$true <- pd$beta[S0][res$coord]
  res$full <- unname(f$beta[match(S0, Ai)])[res$coord]
  list(gate_fail = FALSE, res = res, diag = c(events_inf = sum(d$status[inf]), n_inf = length(inf)))
}

## ---------------------------------------------------------------------
## Plasmode study driver (run directory with app/ processed cohorts and
## app/results from summarize_apps()).  For each cohort, the true support is
## the set of variables selected in at least 95 of the 100 splits of the
## data analysis and the refitted set is the primary-split selected set.
## ---------------------------------------------------------------------
PLASMODE_R <- 500L
PLASMODE_CHUNK <- 25L
PLASMODE_SEED <- 20261026L

prepare_plasmode <- function(dir = ".", cohorts = names(app_cohorts()), stable = 0.95) {
  od <- file.path(dir, "app", "plasmode"); dir.create(od, recursive = TRUE, showWarnings = FALSE)
  P <- list()
  for (co in cohorts) {
    ff <- file.path(dir, "app", "results", paste0(co, c("_selection_frequency.csv", "_primary.csv")))
    if (!all(file.exists(ff)) || !file.exists(file.path(dir, app_cohorts()[[co]]$file))) next
    fr <- utils::read.csv(ff[1]); A <- utils::read.csv(ff[2])$var
    sup <- intersect(A, fr$var[fr$frequency >= stable])
    d <- load_cohort(dir, co)
    pd <- plasmode_design(d$data, d$covs, sup, tau = d$tau)
    pd$X <- NULL                         # rebuilt from the cohort file in each task
    P[[co]] <- list(pd = pd, A = A, support = sup, r = app_cohorts()[[co]]$r)
    message(sprintf("%-16s |S0| = %d, |A| = %d, lambda0 = %.4g (crude %.4g), events %.3f -> %.3f", co, length(sup),
                    length(A), pd$lambda0, pd$crude, pd$ev_obs, pd$ev_sim))
  }
  saveRDS(P, file.path(od, "plasmode_designs.rds")); invisible(P)
}

plasmode_task_table <- function(dir = ".") {
  P <- readRDS(file.path(dir, "app", "plasmode", "plasmode_designs.rds"))
  first <- seq(1L, PLASMODE_R, by = PLASMODE_CHUNK)
  tt <- expand.grid(first = first, cohort = names(P), stringsAsFactors = FALSE)
  data.frame(task = seq_len(nrow(tt)), cohort = tt$cohort, first = tt$first,
             last = pmin(tt$first + PLASMODE_CHUNK - 1L, PLASMODE_R), stringsAsFactors = FALSE)
}

run_plasmode_task <- function(task, dir = ".") {
  P <- readRDS(file.path(dir, "app", "plasmode", "plasmode_designs.rds")); tt <- plasmode_task_table(dir)
  tk <- tt[task, ]; cfg <- P[[tk$cohort]]
  od <- file.path(dir, "app", "plasmode", tk$cohort); dir.create(od, showWarnings = FALSE)
  of <- file.path(od, sprintf("reps_%04d_%04d.rds", tk$first, tk$last))
  if (file.exists(of)) { message("exists: ", of); return(invisible(of)) }
  pd <- cfg$pd; pd$X <- load_cohort(dir, tk$cohort)$X
  ci <- match(tk$cohort, names(P)); res <- list()
  for (k in tk$first:tk$last) {
    t0 <- tic()
    res[[as.character(k)]] <- tryCatch(run_rep_plasmode_known(pd, cfg$A, PLASMODE_SEED + 10000L * ci + k, cfg$r),
                                       error = function(e) list(error = conditionMessage(e)))
    message(sprintf("[%s] rep %d %.1fs", tk$cohort, k, tic() - t0))
  }
  saveRDS(list(cohort = tk$cohort, reps = res), of); invisible(of)
}

summarize_plasmode <- function(dir = ".") {
  P <- readRDS(file.path(dir, "app", "plasmode", "plasmode_designs.rds"))
  od <- file.path(dir, "app", "results"); co_l <- list(); av_l <- list()
  for (co in names(P)) {
    fs <- list.files(file.path(dir, "app", "plasmode", co), "^reps_.*\\.rds$", full.names = TRUE)
    if (!length(fs)) next
    reps <- unlist(lapply(fs, function(f) readRDS(f)$reps), recursive = FALSE)
    cs <- summ_coord(reps, co); co_l[[co]] <- cs; av_l[[co]] <- summ_avg(cs)
  }
  cs <- do.call(rbind, co_l); av <- do.call(rbind, av_l)
  info <- do.call(rbind, lapply(names(P), function(co) with(P[[co]], data.frame(cohort = co, n_support = length(support),
    n_A = length(A), lambda0 = pd$lambda0, crude = pd$crude, events_obs = pd$ev_obs, events_sim = pd$ev_sim,
    support = paste(support, collapse = " ")))))
  utils::write.csv(cs, file.path(od, "plasmode_coordinate_level.csv"), row.names = FALSE)
  utils::write.csv(av, file.path(od, "plasmode_averaged.csv"), row.names = FALSE)
  utils::write.csv(info, file.path(od, "plasmode_designs.csv"), row.names = FALSE)
  lab <- c(msk_met = "MSK-MET clinical", msk_met_genomic = "MSK-MET genomic", nhis_lmf = "NHIS-LMF")
  fig_plasmode2(av, file.path(od, "fig_plasmode.pdf"), lab[intersect(names(lab), unique(av$design))])
  invisible(list(coord = cs, avg = av, info = info))
}
