## =====================================================================
## File 05_simrep.R -- one replicate of each simulation study.
## =====================================================================

## Methods evaluated at the refit stage.  name = list(rule, score, C)
ipw_methods_default <- function() list(
  UIPW = list(rule = "uniform", score = "q"),
  SIPW = list(rule = "scorenorm", score = "q"),
  ZIPW = list(rule = "zuo", score = "q"),
  IIPW = list(rule = "inference", score = "q"),
  CIPW = list(rule = "inference", score = "q", contrast = 1L),   # targets coordinate 1 only
  ## Poisson (Bernoulli) designs, events included with certainty when #events < r
  PUIPW = list(design = "poisson", rule = "uniform"),
  PIIPW = list(design = "poisson", rule = "inference"))

## one Poisson method at all budgets
.run_poisson <- function(mn, mm, st, pil, r_vec, Up, A, S0, t_base, C = NULL) {
  out <- list(); t0 <- tic(); s <- poisson_scores(st, pil, mm$rule, C); t_s <- tic() - t0
  for (r in r_vec) {
    t0 <- tic(); pr <- poisson_inclusion(s, st$ev, r)
    g <- ah_fit_poisson(st, pr$p, Up)
    out[[length(out) + 1]] <- .rec(mn, r, g$beta, g$se, A, S0,
      list(gamma_frac = g$gamma_frac, ess = g$size, max_w = g$max_w, time = t_base + t_s + tic() - t0))
    out[[length(out)]]$event_complete <- pr$event_complete
    ## conditional (subsampling-only) variance of the existing subsampling
    ## literature, which targets the full-sample estimator
    if (mn == "PIIPW") {
      out[[length(out) + 1]] <- .rec("PIIPW_GammaOnly", r, g$beta, sqrt(pmax(diag(g$Sigma_gamma_only), 0)), A, S0)
      out[[length(out)]]$event_complete <- pr$event_complete
    }
    ## Kulich-Lin case-cohort estimator on the same event-complete subsample
    if (mn == "PUIPW") {
      if (pr$event_complete) {
        t0 <- tic(); k <- ah_fit_casecohort(st, pr$p, Up)
        out[[length(out) + 1]] <- .rec("CaseCohort", r, k$beta, k$se, A, S0,
          list(gamma_frac = k$gamma_frac, ess = k$size, max_w = k$max_w, time = t_base + tic() - t0))
      } else out[[length(out) + 1]] <- .rec("CaseCohort", r, numeric(0), numeric(0), integer(0), S0)
      out[[length(out)]]$event_complete <- pr$event_complete
    }
  }
  out
}

.rec <- function(method, r, beta, se, A, S0, extra = list()) {
  pos <- match(S0, A)
  data.frame(method = method, r = r, coord = seq_along(S0),
             est = if (length(beta)) unname(beta[pos]) else NA_real_,
             se = if (length(se)) unname(se[pos]) else NA_real_,
             retained = !is.na(pos),
             gamma_frac = extra$gamma_frac %||% NA_real_, ess = extra$ess %||% NA_real_,
             max_w = extra$max_w %||% NA_real_, time = extra$time %||% NA_real_,
             event_complete = NA, stringsAsFactors = FALSE)
}

## Study 1: support-known (A = S0 plus `n_noise` noise covariates) refit
## stage.  Pilot and inference samples are independent, each of size m.
run_rep_known <- function(des, m, r_vec, seed, n_noise = 4, methods = ipw_methods_default(),
                          resid = TRUE, target = NULL) {
  set.seed(seed)
  S0 <- des$S0; A <- c(S0, setdiff(seq_len(des$p), S0)[seq_len(n_noise)])
  dp <- gen_data(des, n = m); di <- gen_data(des, n = m)
  if (mean(di$time >= des$tau) < 0.05) return(list(gate_fail = TRUE))
  t0 <- tic(); pil <- pilot_fit(dp$time, dp$status, dp$X[, A, drop = FALSE], des$tau); t_pil <- tic() - t0
  t0 <- tic(); st <- ah_setup(di$time, di$status, di$X[, A, drop = FALSE], des$tau); t_st <- tic() - t0
  out <- list()
  t0 <- tic(); f <- ah_fit_full(st)
  out[[1]] <- .rec("SplitFull", NA, f$beta, f$se, A, S0, list(time = t_st + tic() - t0))
  U <- stats::runif(max(r_vec)); Up <- stats::runif(st$n)
  for (mn in names(methods)) {
    mm <- methods[[mn]]
    C <- if (!is.null(mm$contrast)) matrix(as.numeric(seq_along(A) == mm$contrast), 1) else NULL
    if (identical(mm$design, "poisson")) {
      out <- c(out, .run_poisson(mn, mm, st, pil, r_vec, Up, A, S0, t_pil + t_st, C)); next
    }
    t0 <- tic(); pi <- ipw_probs(st, pil, mm$rule, score = mm$score, C = C); t_pi <- tic() - t0
    for (r in r_vec) {
      t0 <- tic(); g <- ah_fit_ipw(st, pi, draw_ipw(pi, U[seq_len(r)]), score = mm$score)
      out[[length(out) + 1]] <- .rec(mn, r, g$beta, g$se, A, S0,
        list(gamma_frac = g$gamma_frac, ess = g$ess, max_w = g$max_w, time = t_pil + t_st + t_pi + tic() - t0))
      if (mn == "IIPW")   # variance ablation: same estimate, Omega-only SE
        out[[length(out) + 1]] <- .rec("IIPW_OmegaOnly", r, g$beta, sqrt(diag(g$Sigma_omega_only)), A, S0)
      if (mn %in% c("IIPW", "ZIPW"))   # conditional variance of the subsampling literature
        out[[length(out) + 1]] <- .rec(paste0(mn, "_GammaOnly"), r, g$beta, sqrt(diag(g$Sigma_gamma_only)), A, S0)
    }
  }
  if (resid) {             # exploratory residual-score IPW (not in the manuscript)
    pi <- ipw_probs(st, pil, "inference", score = "xi")
    for (r in r_vec) { g <- ah_fit_ipw(st, pi, draw_ipw(pi, U[seq_len(r)]), score = "xi")
      out[[length(out) + 1]] <- .rec("RIPW", r, g$beta, g$se, A, S0, list(gamma_frac = g$gamma_frac)) }
  }
  res <- do.call(rbind, out); res$true <- (target %||% des$beta0[S0])[res$coord]
  res$full <- f$beta[match(S0, A)][res$coord]        # full-inference-sample estimate
  list(gate_fail = FALSE, res = res,
       diag = c(cens = mean(di$status == 0), event_tau = mean(st$ev), risk_tau = mean(di$time >= des$tau)))
}

## Studies 2-3: full pipeline (split, screen, pilot, refits).
## screen: "gselo" | "marginal" (direct top-K) | "marginal_gselo"
run_rep_pipeline <- function(des, seed, r_vec = des$r, methods = ipw_methods_default()[c("UIPW", "IIPW")],
                             screen = "gselo", K = 25, naive = FALSE, oracle = TRUE, target = NULL) {
  set.seed(seed)
  d <- gen_data(des)
  n <- length(d$time)
  if (mean(d$time >= des$tau) < 0.05) return(list(gate_fail = TRUE))
  scr_idx <- sort(sample.int(n, floor(n / 2))); inf_idx <- setdiff(seq_len(n), scr_idx)
  S0 <- des$S0; tgt <- target %||% des$beta0
  ## screening
  sc <- switch(screen,
    gselo = screen_gselo(d$time[scr_idx], d$status[scr_idx], d$X[scr_idx, , drop = FALSE], des$tau),
    marginal = screen_marginal(d$time[scr_idx], d$status[scr_idx], d$X[scr_idx, , drop = FALSE], des$tau, K),
    marginal_gselo = screen_marginal(d$time[scr_idx], d$status[scr_idx], d$X[scr_idx, , drop = FALSE],
                                     des$tau, K, then_gselo = TRUE))
  A <- sc$A
  sdiag <- c(size = length(A), sure = all(S0 %in% A), exact = setequal(A, S0),
             theta = sc$theta %||% NA, eta = sc$eta %||% NA,
             b_low = sc$boundary_low %||% NA, b_high = sc$boundary_high %||% NA,
             kkt = sc$kkt %||% NA, conv = sc$all_converged %||% NA, t_screen = sc$time,
             cand_sure = if (!is.null(sc$candidates)) all(S0 %in% sc$candidates) else NA)
  out <- list()
  if (length(A)) {
    Xs <- d$X[scr_idx, A, drop = FALSE]; Xi <- d$X[inf_idx, A, drop = FALSE]
    t0 <- tic(); pil <- pilot_fit(d$time[scr_idx], d$status[scr_idx], Xs, des$tau); t_pil <- tic() - t0
    t0 <- tic(); st <- ah_setup(d$time[inf_idx], d$status[inf_idx], Xi, des$tau); t_st <- tic() - t0
    f <- ah_fit_full(st)
    out[[1]] <- .rec("SplitFull", NA, f$beta, f$se, A, S0, list(time = sc$time + t_st))
    U <- stats::runif(max(r_vec)); Up <- stats::runif(st$n)
    for (mn in names(methods)) {
      mm <- methods[[mn]]
      C <- if (!is.null(mm$contrast)) matrix(as.numeric(A == S0[mm$contrast]), 1) else NULL
      if (identical(mm$design, "poisson")) {
        out <- c(out, .run_poisson(mn, mm, st, pil, r_vec, Up, A, S0, sc$time + t_pil + t_st, C)); next
      }
      t0 <- tic(); pi <- ipw_probs(st, pil, mm$rule, score = mm$score, C = C); t_pi <- tic() - t0
      for (r in r_vec) {
        t0 <- tic(); g <- ah_fit_ipw(st, pi, draw_ipw(pi, U[seq_len(r)]), score = mm$score)
        out[[length(out) + 1]] <- .rec(mn, r, g$beta, g$se, A, S0,
          list(gamma_frac = g$gamma_frac, ess = g$ess, max_w = g$max_w,
               time = sc$time + t_pil + t_st + t_pi + tic() - t0))
        if (mn %in% c("IIPW", "ZIPW"))
          out[[length(out) + 1]] <- .rec(paste0(mn, "_GammaOnly"), r, g$beta, sqrt(diag(g$Sigma_gamma_only)), A, S0)
      }
    }
  } else {
    out[[1]] <- .rec("SplitFull", NA, numeric(0), numeric(0), integer(0), S0)
    for (mn in names(methods)) for (r in r_vec) {
      out[[length(out) + 1]] <- .rec(mn, r, numeric(0), numeric(0), integer(0), S0)
      if (mn == "PUIPW") out[[length(out) + 1]] <- .rec("CaseCohort", r, numeric(0), numeric(0), integer(0), S0)
    }
  }
  if (oracle) {     # support-known refit on the inference split (no screening)
    st0 <- ah_setup(d$time[inf_idx], d$status[inf_idx], d$X[inf_idx, S0, drop = FALSE], des$tau)
    f0 <- ah_fit_full(st0)
    out[[length(out) + 1]] <- .rec("OracleFull", NA, f0$beta, f0$se, S0, S0)
  }
  if (naive) {      # screen and refit on the same full data (no splitting)
    t0 <- tic()
    scn <- screen_gselo(d$time, d$status, d$X, des$tau)
    if (length(scn$A)) {
      stn <- ah_setup(d$time, d$status, d$X[, scn$A, drop = FALSE], des$tau)
      fn <- ah_fit_full(stn)
      out[[length(out) + 1]] <- .rec("Naive", NA, fn$beta, fn$se, scn$A, S0, list(time = tic() - t0))
    } else out[[length(out) + 1]] <- .rec("Naive", NA, numeric(0), numeric(0), integer(0), S0)
  }
  res <- do.call(rbind, out); res$true <- tgt[S0][res$coord]
  fb <- rep(NA_real_, length(S0))
  if (length(A)) { pos <- match(S0, A); fb[!is.na(pos)] <- f$beta[pos[!is.na(pos)]] }
  res$full <- fb[res$coord]
  list(gate_fail = FALSE, res = res, screen = sdiag,
       diag = c(cens = mean(d$status == 0), risk_tau = mean(d$time >= des$tau)))
}

## Study 3 (p >= n_scr): marginal pre-screen candidate retention for several
## K and the optional candidate-restricted GSELO reduction
run_rep_marginal <- function(des, seed, Ks = c(25, 50)) {
  set.seed(seed)
  d <- gen_data(des); n <- length(d$time)
  scr <- sort(sample.int(n, floor(n / 2)))
  bv <- ah_bV(d$time[scr], d$status[scr], d$X[scr, , drop = FALSE], des$tau)
  u <- abs(bv$b) / sqrt(pmax(bv$diagV, 1e-300))
  o <- order(-u, seq_along(u)); S0 <- des$S0
  out <- NULL
  for (K in Ks) {
    cand <- sort(o[seq_len(K)]); csure <- all(S0 %in% cand)
    g <- tryCatch(screen_gselo(d$time[scr], d$status[scr], d$X[scr, , drop = FALSE], des$tau,
                               candidates = cand), error = function(e) NULL)
    out <- rbind(out, data.frame(K = K, cand_sure = csure,
      gselo_sure = if (is.null(g)) NA else all(S0 %in% g$A), gselo_fail = is.null(g),
      gselo_size = if (is.null(g)) NA else length(g$A), risk_tau = mean(d$time >= des$tau)))
  }
  list(gate_fail = FALSE, res = out)
}
