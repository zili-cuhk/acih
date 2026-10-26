## =====================================================================
## File 06_realdata.R -- the data-analysis protocol of Section 5 and Web
## Appendix C, for any processed cohort (time, status, covariate matrix).
##  * prespecified primary 50/50 split (seed 20260701) and 99 further splits;
##  * primary screen: GSELO continuation path + joint EBIC (as in the
##    simulations); fixed-budget screens (top-K |coefficients| of the
##    unpenalised multivariable fit) with K = 5, 10, 15 as sensitivity;
##  * split-full refit with martingale sandwich (primary table);
##  * with-replacement IPW (uniform, score norm, inference-oriented with
##    C = I), Poisson designs (event-complete when #events < r; uniform or
##    inference-oriented censored probabilities) and the Kulich-Lin
##    case-cohort estimator on the PUIPW subsample, at several budgets;
##  * fixed-support survey-weighted sensitivity with stratified PSU
##    linearisation (NHIS-LMF); Kaplan-Meier curves by score tertile.
## One split of one cohort is one task (LSF array index), so the 100 splits
## run in parallel; `summarize_apps()` collects them.
## =====================================================================

primary_split <- function(n, seed) {
  set.seed(seed); s <- sort(sample.int(n, floor(n / 2))); list(scr = s, inf = setdiff(seq_len(n), s))
}

APP_SEED <- 20260701L
APP_R <- c(1000, 2000, 5000, 10000)

## cohorts analysed in Section 5: processed file (relative to the run
## directory), horizon tau and its unit
app_cohorts <- function() list(
  msk_met         = list(file = "app/msk_met_clinical.rds", tau = 60, unit = "month", label = "MSK-MET clinical",
                         r = c(1000, 2000, 7500, 10000), cap = 40),
  msk_met_genomic = list(file = "app/msk_met_genomic.rds",  tau = 60, unit = "month", label = "MSK-MET genomic",
                         r = c(1000, 2000, 7500, 10000), cap = 40),
  nhis_lmf        = list(file = "app/nhis_lmf.rds",         tau = 3,  unit = "year",  label = "NHIS-LMF",
                         r = c(1000, 2000, 5000, 10000), cap = 40))
## size cap of the GSELO screen: 40 in the applications (20 in the
## simulations); `cap = 20` is reported as a sensitivity analysis
## budgets: about 5,000 of the 12,281 MSK-MET inference records are deaths
## within 60 months, so the event-complete design is used at r = 7500 and
## 10000 (r = 5000 would leave only a few censored records); about 1,270 of
## the 41,838 NHIS-LMF inference records are deaths within three years.

## ---------------------------------------------------------------------
## one split: screen on D_scr, refits on D_inf
## ---------------------------------------------------------------------
analyze_split <- function(time, status, X, tau, split, screen = c("gselo", "topk"), budget = 10, cap = 20,
                          r_vec = APP_R, seed_draw = 1, subsample = TRUE) {
  screen <- match.arg(screen)
  s <- split$scr; i <- split$inf; vn <- colnames(X) %||% paste0("X", seq_len(ncol(X)))
  t0 <- tic()
  sc <- if (screen == "topk") screen_topk_unpen(time[s], status[s], X[s, , drop = FALSE], tau, budget)
        else screen_gselo(time[s], status[s], X[s, , drop = FALSE], tau, cap = cap)
  t_scr <- tic() - t0
  A <- sort(sc$A)
  info <- data.frame(screen = screen, size = length(A), t_screen = t_scr,
                     theta = sc$theta %||% NA_real_, eta = sc$eta %||% NA_real_,
                     boundary_low = sc$boundary_low %||% NA, boundary_high = sc$boundary_high %||% NA,
                     converged = sc$all_converged %||% NA,
                     events_scr = sum(status[s] == 1 & time[s] <= tau), events_inf = sum(status[i] == 1 & time[i] <= tau),
                     n_inf = length(i))
  if (!length(A)) return(list(A = character(0), info = info, table = NULL))
  rec <- function(method, r, g, tt, extra = list())
    data.frame(method = method, r = r, var = vn[A], est = g$beta, se = g$se, time = tt,
               gamma_frac = extra$gamma_frac %||% NA_real_, size = extra$size %||% NA_real_,
               event_complete = extra$event_complete %||% NA, stringsAsFactors = FALSE)
  t0 <- tic(); st <- ah_setup(time[i], status[i], X[i, A, drop = FALSE], tau); t_st <- tic() - t0
  t0 <- tic(); full <- ah_fit_full(st); t_full <- tic() - t0
  out <- list(rec("SplitFull", NA_real_, full, t_scr + t_st + t_full, list(size = st$n)))
  ## diagnostic: fraction of inference-split records whose fitted hazard is
  ## negative at the average Breslow baseline over [0, tau]
  lp <- drop(X[i, A, drop = FALSE] %*% full$beta)
  lam0 <- sum(st$ob$dk / st$ob$S0) / tau - mean(st$Xbar %*% full$beta)
  info$neg_hazard <- mean(lam0 + lp < 0)
  if (subsample && length(r_vec)) {
    t0 <- tic(); pil <- pilot_fit(time[s], status[s], X[s, A, drop = FALSE], tau); t_pil <- tic() - t0
    t_base <- t_scr + t_st + t_pil
    set.seed(seed_draw); U <- stats::runif(max(r_vec)); Up <- stats::runif(st$n)
    for (rule in c(UIPW = "uniform", SIPW = "scorenorm", IIPW = "inference")) {
      mn <- names(which(c(UIPW = "uniform", SIPW = "scorenorm", IIPW = "inference") == rule))
      t0 <- tic(); pi <- ipw_probs(st, pil, rule); t_pi <- tic() - t0
      for (r in r_vec) {
        t0 <- tic(); g <- ah_fit_ipw(st, pi, draw_ipw(pi, U[seq_len(r)]))
        out[[length(out) + 1]] <- rec(mn, r, g, t_base + t_pi + tic() - t0, list(gamma_frac = g$gamma_frac, size = r))
      }
    }
    for (rule in c(PUIPW = "uniform", PIIPW = "inference")) {
      mn <- names(which(c(PUIPW = "uniform", PIIPW = "inference") == rule))
      t0 <- tic(); sp <- poisson_scores(st, pil, rule); t_ps <- tic() - t0
      for (r in r_vec) {
        t0 <- tic(); pr <- poisson_inclusion(sp, st$ev, r); g <- ah_fit_poisson(st, pr$p, Up)
        out[[length(out) + 1]] <- rec(mn, r, g, t_base + t_ps + tic() - t0,
                                      list(gamma_frac = g$gamma_frac, size = g$size, event_complete = pr$event_complete))
        if (mn == "PUIPW" && pr$event_complete) {     # Kulich-Lin on the same subsample
          t0 <- tic(); k <- ah_fit_casecohort(st, pr$p, Up)
          out[[length(out) + 1]] <- rec("CaseCohort", r, k, t_scr + tic() - t0,
                                        list(gamma_frac = k$gamma_frac, size = k$size, event_complete = TRUE))
        }
      }
    }
  }
  list(A = vn[A], info = info, table = do.call(rbind, out))
}

## Survey-weighted fixed-support fit with stratified PSU linearisation
ah_fit_survey <- function(time, status, X, tau, weight, strata, psu) {
  w <- weight / mean(weight)
  st <- ah_setup(time, status, X, tau, w = w)
  beta <- solve_ridge(st$Vfull, st$bfull)
  xi <- xi_scores(st, beta) * st$w                  # weighted subject scores (sorted order)
  o <- st$ob$o; h <- strata[o]; g <- psu[o]
  Q <- rowsum(xi, paste(h, g, sep = "\r"))
  hs <- sub("\r.*$", "", rownames(Q))
  meat <- matrix(0, st$d, st$d); m_h <- integer(0)
  for (hh in unique(hs)) {
    Qh <- Q[hs == hh, , drop = FALSE]; mh <- nrow(Qh); m_h <- c(m_h, mh); if (mh < 2) next
    Qc <- sweep(Qh, 2, colMeans(Qh)); meat <- meat + mh / (mh - 1) * crossprod(Qc)
  }
  meat <- meat / st$W^2
  Hi <- solve(st$Vfull)
  V <- Hi %*% meat %*% t(Hi)
  list(beta = beta, se = sqrt(diag(V)), n_strata = length(unique(hs)), n_psu = nrow(Q), min_psu = min(m_h))
}

## Kaplan-Meier (base R) for descriptive curves
km <- function(time, status) {
  o <- order(time); t <- time[o]; s <- status[o]
  ut <- sort(unique(t[s == 1]))
  atrisk <- length(t) - findInterval(ut, t, left.open = TRUE)
  dd <- as.numeric(table(factor(t[s == 1], levels = ut)))
  list(time = ut, surv = cumprod(1 - dd / atrisk))
}

## ---------------------------------------------------------------------
## task table and task runner (one cohort x one split per task)
## ---------------------------------------------------------------------
app_task_table <- function(dir = ".", n_splits = 100, cohorts = app_cohorts()) {
  have <- vapply(cohorts, function(x) file.exists(file.path(dir, x$file)), TRUE)
  if (!any(have)) stop("no processed cohort found; run the app/*_prepare.R scripts first")
  tt <- expand.grid(split = seq_len(n_splits), cohort = names(cohorts)[have], stringsAsFactors = FALSE)
  data.frame(task = seq_len(nrow(tt)), cohort = tt$cohort, split = tt$split, stringsAsFactors = FALSE)
}

load_cohort <- function(dir, cohort) {
  cf <- app_cohorts()[[cohort]]; d <- readRDS(file.path(dir, cf$file))
  X <- as.matrix(d$data[, d$covs, drop = FALSE]); colnames(X) <- d$covs
  c(d, list(X = X, tau = cf$tau, cfg = cf))
}

run_app_task <- function(task, dir = ".", n_splits = 100, r_vec = NULL) {
  tt <- app_task_table(dir, n_splits); if (task < 1 || task > nrow(tt)) stop("no such task: ", task)
  tk <- tt[task, ]
  od <- file.path(dir, "app", "out", tk$cohort); dir.create(od, recursive = TRUE, showWarnings = FALSE)
  of <- file.path(od, sprintf("split_%03d.rds", tk$split))
  if (file.exists(of)) { message("exists: ", of); return(invisible(of)) }
  co <- load_cohort(dir, tk$cohort); dat <- co$data; n <- nrow(co$X)
  if (is.null(r_vec)) r_vec <- co$cfg$r %||% APP_R
  sp <- primary_split(n, if (tk$split == 1) APP_SEED else APP_SEED + tk$split)
  t0 <- tic()
  res <- list(primary = analyze_split(dat$time, dat$status, co$X, co$tau, sp, "gselo", cap = co$cfg$cap %||% 20,
                                              r_vec = r_vec, seed_draw = tk$split))
  if (tk$split == 1) {
    ## fixed-budget sensitivity on the primary split (split-full only)
    res$budget <- lapply(c(5, 10, 15), function(K)
      analyze_split(dat$time, dat$status, co$X, co$tau, sp, "topk", budget = min(K, ncol(co$X)), subsample = FALSE))
    names(res$budget) <- c("K5", "K10", "K15")
    ## fixed-support survey sensitivity on the full processed cohort
    if (!is.null(co$survey) && length(res$primary$A)) {
      XA <- co$X[, res$primary$A, drop = FALSE]
      u <- ah_fit_survey(dat$time, dat$status, XA, co$tau, rep(1, n), co$survey$strata, co$survey$psu)
      w <- ah_fit_survey(dat$time, dat$status, XA, co$tau, co$survey$weight, co$survey$strata, co$survey$psu)
      res$survey <- data.frame(var = res$primary$A, U_est = u$beta, U_se = u$se, W_est = w$beta, W_se = w$se)
      res$survey_design <- c(n_strata = w$n_strata, n_psu = w$n_psu, min_psu = w$min_psu)
    }
  }
  res$cohort <- tk$cohort; res$split <- tk$split; res$seconds <- tic() - t0; res$R_version <- R.version.string
  saveRDS(res, of)
  message(sprintf("[%s split %d] |A| = %d, %.1f s", tk$cohort, tk$split, length(res$primary$A), res$seconds))
  invisible(of)
}

## ---------------------------------------------------------------------
## summaries, tables and figures
## ---------------------------------------------------------------------
summarize_apps <- function(dir = ".", cohorts = names(app_cohorts())) {
  od <- file.path(dir, "app", "results"); dir.create(od, showWarnings = FALSE, recursive = TRUE)
  S <- list()
  for (co in cohorts) {
    fs <- sort(list.files(file.path(dir, "app", "out", co), "^split_.*\\.rds$", full.names = TRUE))
    if (!length(fs)) next
    L <- lapply(fs, readRDS); K <- length(L)
    pr <- L[[which(vapply(L, function(x) x$split, 1) == 1)]]
    info <- do.call(rbind, lapply(L, function(x) cbind(split = x$split, x$primary$info)))
    tab <- do.call(rbind, lapply(L, function(x) if (!is.null(x$primary$table)) cbind(split = x$split, x$primary$table)))
    sf <- tab[tab$method == "SplitFull", ]
    allv <- colnames(load_cohort(dir, co)$X)
    freq <- table(factor(sf$var, levels = allv)) / K
    ## primary table with selection frequency and sign agreement over splits
    p1 <- pr$primary$table[pr$primary$table$method == "SplitFull", c("var", "est", "se")]
    p1$lo <- p1$est - 1.96 * p1$se; p1$hi <- p1$est + 1.96 * p1$se
    p1$frequency <- as.numeric(freq[p1$var])
    p1$sign_agree <- vapply(seq_len(nrow(p1)), function(j) {
      x <- sf$est[sf$var == p1$var[j]]; mean(sign(x) == sign(p1$est[j])) }, 1)
    ## efficiency: per split, mean SE over the selected coordinates and its
    ## ratio to SplitFull; then averaged over splits
    tab$key <- paste(tab$split, tab$var)
    sfse <- stats::setNames(sf$se, paste(sf$split, sf$var))
    tab$se_rel <- tab$se / sfse[tab$key]
    agg <- stats::aggregate(cbind(se, se_rel, time, gamma_frac, size) ~ method + r + split,
                            data = transform_r(tab), FUN = mean, na.action = stats::na.pass)
    eff <- stats::aggregate(cbind(se, se_rel, time, gamma_frac, size) ~ method + r, data = agg, FUN = mean,
                            na.action = stats::na.pass)
    nspl <- stats::aggregate(split ~ method + r, data = agg, FUN = length); names(nspl)[3] <- "n_splits"
    eff <- merge(eff, nspl); eff$r[eff$r == 0] <- NA
    eff <- eff[order(match(eff$method, METHODS_APP), eff$r), ]
    sel <- data.frame(cohort = co, splits = K, size_mean = mean(info$size), size_min = min(info$size),
                      size_max = max(info$size), cap = mean(info$size >= (app_cohorts()[[co]]$cap %||% 20)), converged = mean(info$converged),
                      boundary = mean(info$boundary_low | info$boundary_high), t_screen = mean(info$t_screen),
                      events_inf = mean(info$events_inf), n_inf = mean(info$n_inf),
                      neg_hazard = mean(info$neg_hazard, na.rm = TRUE))
    bud <- NULL
    if (!is.null(pr$budget)) {
      vv <- unique(unlist(lapply(pr$budget, `[[`, "A")))
      bud <- data.frame(var = vv)
      for (nm in names(pr$budget)) { t1 <- pr$budget[[nm]]$table; t1 <- t1[t1$method == "SplitFull", ]
        bud[[nm]] <- t1$est[match(vv, t1$var)]; bud[[paste0(nm, "_se")]] <- t1$se[match(vv, t1$var)] }
      g <- p1$est[match(vv, p1$var)]; bud$GSELO <- g
    }
    utils::write.csv(p1, file.path(od, paste0(co, "_primary.csv")), row.names = FALSE)
    utils::write.csv(eff, file.path(od, paste0(co, "_efficiency.csv")), row.names = FALSE)
    utils::write.csv(data.frame(var = allv, frequency = as.numeric(freq)), file.path(od, paste0(co, "_selection_frequency.csv")), row.names = FALSE)
    utils::write.csv(sel, file.path(od, paste0(co, "_screening.csv")), row.names = FALSE)
    if (!is.null(bud)) utils::write.csv(bud, file.path(od, paste0(co, "_budget_sensitivity.csv")), row.names = FALSE)
    if (!is.null(pr$survey)) utils::write.csv(cbind(pr$survey, t(pr$survey_design)), file.path(od, paste0(co, "_survey_sensitivity.csv")), row.names = FALSE)
    S[[co]] <- list(primary = p1, efficiency = eff, screening = sel, budget = bud, survey = pr$survey,
                    survey_design = pr$survey_design, frequency = freq, splitfull = sf, info = info)
  }
  saveRDS(S, file.path(od, "app_summary.rds"))
  if (length(S)) {
    lab <- c(msk_met = "MSK-MET clinical (p = 20)", msk_met_genomic = "MSK-MET genomic (p = 223)", nhis_lmf = "NHIS-LMF (3% mortality)")
    fig_app_efficiency2(lapply(S, `[[`, "efficiency"), file.path(od, "fig_real_data_efficiency.pdf"), lab[names(S)])
    km <- app_km_data(dir, S)
    utils::write.csv(km$curves, file.path(od, "km_curves.csv"), row.names = FALSE)
    utils::write.csv(km$risk, file.path(od, "km_risk.csv"), row.names = FALSE)
    fig_km2(km$curves, km$risk, file.path(od, "fig_real_data_km.pdf"), lab[names(S)],
            vapply(names(S), function(x) paste0(app_cohorts()[[x]]$unit, "s"), ""))
  }
  ## size-cap sensitivity of the genomic analysis (Web Appendix C)
  if ("msk_met_genomic" %in% names(S)) app_cap_sensitivity(dir, "msk_met_genomic", cap = 20)
  invisible(S)
}
METHODS_APP <- c("SplitFull", "UIPW", "SIPW", "IIPW", "PUIPW", "PIIPW", "CaseCohort")
transform_r <- function(tab) { tab$r[is.na(tab$r)] <- 0; tab }

## SE relative to the split-full refit against the budget r (log scale)
fig_app_efficiency <- function(S, file, methods = c("UIPW", "SIPW", "IIPW", "PUIPW", "PIIPW", "CaseCohort")) {
  col <- c(UIPW = "#1b9e77", SIPW = "#7570b3", IIPW = "#d95f02", PUIPW = "#6a3d9a", PIIPW = "#e7298a", CaseCohort = "#666666")
  grDevices::pdf(file, width = 3.4 * length(S) + 0.4, height = 3.6)
  op <- graphics::par(mfrow = c(1, length(S)), mar = c(4, 4.2, 2, 0.6), oma = c(2.2, 0, 0, 0)); on.exit({ graphics::par(op); grDevices::dev.off() })
  lab <- vapply(names(S), function(x) app_cohorts()[[x]]$label, "")
  for (k in seq_along(S)) {
    e <- S[[k]]$efficiency; e <- e[e$method %in% methods, ]
    rr <- sort(unique(e$r)); yl <- range(c(1, e$se_rel), na.rm = TRUE) * c(0.9, 1.1)
    graphics::plot(NA, xlim = range(rr), ylim = yl, log = "xy", xaxt = "n", xlab = "Subsample size r",
                   ylab = "Mean SE / split-full SE", main = sprintf("(%s) %s", LETTERS[k], lab[k]), cex.main = 0.95, font.main = 1)
    graphics::axis(1, at = rr, labels = format(rr, big.mark = ","), cex.axis = 0.8)
    graphics::abline(h = 1, lty = 2)
    for (m in methods) { x <- e[e$method == m, ]; if (!nrow(x)) next
      graphics::lines(x$r, x$se_rel, col = col[m], lwd = 2); graphics::points(x$r, x$se_rel, col = col[m], pch = 19, cex = 0.8) }
  }
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  graphics::legend("bottom", legend = methods, col = col[methods], lwd = 2, pch = 19, horiz = TRUE, bty = "n", cex = 0.8, xpd = TRUE)
}

## Kaplan-Meier curves by tertile of the primary split-full linear
## predictor (same-data description)
fig_app_km <- function(dir, S, file) {
  grDevices::pdf(file, width = 3.4 * length(S) + 0.4, height = 3.6)
  op <- graphics::par(mfrow = c(1, length(S)), mar = c(4, 4.2, 2, 0.6)); on.exit({ graphics::par(op); grDevices::dev.off() })
  km_tab <- NULL
  for (k in seq_along(S)) {
    co <- names(S)[k]; d <- load_cohort(dir, co); p1 <- S[[co]]$primary
    lp <- drop(d$X[, p1$var, drop = FALSE] %*% p1$est)
    grp <- cut(lp, stats::quantile(lp, c(0, 1 / 3, 2 / 3, 1)), include.lowest = TRUE, labels = c("Low", "Middle", "High"))
    tt <- pmin(d$data$time, d$tau); ss <- d$data$status * (d$data$time <= d$tau)
    cols <- c("#1b9e77", "#7570b3", "#d95f02"); ymin <- 1
    fits <- lapply(levels(grp), function(g) km(tt[grp == g], ss[grp == g]))
    ymin <- min(vapply(fits, function(f) min(c(1, f$surv)), 1))
    graphics::plot(NA, xlim = c(0, d$tau), ylim = c(max(0, ymin - 0.02), 1), xlab = sprintf("Time (%ss)", d$cfg$unit),
                   ylab = "Survival probability", main = sprintf("(%s) %s", LETTERS[k], d$cfg$label), cex.main = 0.95, font.main = 1)
    for (j in 1:3) graphics::lines(c(0, fits[[j]]$time), c(1, fits[[j]]$surv), type = "s", col = cols[j], lwd = 2)
    graphics::legend("bottomleft", legend = paste(levels(grp), "score"), col = cols, lwd = 2, bty = "n", cex = 0.8)
    at <- d$tau * 0.6
    km_tab <- rbind(km_tab, data.frame(cohort = co, group = levels(grp), at = at,
      surv = vapply(fits, function(f) { j <- findInterval(at, f$time); if (j == 0) 1 else f$surv[j] }, 1)))
  }
  utils::write.csv(km_tab, sub("\\.pdf$", ".csv", file), row.names = FALSE)
}

## Kaplan-Meier curves by tertile of the primary-split linear predictor and
## numbers at risk at five time points (input of fig_km2)
app_km_data <- function(dir, S) {
  cu <- list(); ri <- list()
  for (co in names(S)) {
    d <- load_cohort(dir, co); p1 <- S[[co]]$primary
    lp <- drop(d$X[, p1$var, drop = FALSE] %*% p1$est)
    grp <- cut(lp, stats::quantile(lp, c(0, 1 / 3, 2 / 3, 1)), include.lowest = TRUE, labels = c("Low", "Middle", "High"))
    tt <- pmin(d$data$time, d$tau); ss <- d$data$status * (d$data$time <= d$tau)
    at <- pretty(c(0, d$tau), 4); at <- at[at <= d$tau]
    for (g in levels(grp)) {
      f <- km(tt[grp == g], ss[grp == g])
      cu[[length(cu) + 1]] <- data.frame(cohort = co, group = g, time = f$time, surv = f$surv)
      ri[[length(ri) + 1]] <- data.frame(cohort = co, group = g, at = at, n = vapply(at, function(a) sum(tt[grp == g] >= a), 1))
    }
  }
  list(curves = do.call(rbind, cu), risk = do.call(rbind, ri))
}

## ---------------------------------------------------------------------
## Size-cap sensitivity: rerun the GSELO screen of all 100 splits with a
## smaller size cap (20, the simulation cap) and the split-full refit of the
## primary split; writes <cohort>_cap<cap>_frequency.csv (selection
## frequency over splits) and <cohort>_cap<cap>_primary.csv to app/results/
## and returns the proportion of splits whose selected set reached the cap.
## ---------------------------------------------------------------------
app_cap_sensitivity <- function(dir = ".", cohort = "msk_met_genomic", cap = 20, n_split = 100) {
  co <- load_cohort(dir, cohort); n <- nrow(co$X); od <- file.path(dir, "app", "results")
  dir.create(od, showWarnings = FALSE, recursive = TRUE)
  sel <- vector("list", n_split); sz <- numeric(n_split)
  for (k in seq_len(n_split)) {
    sp <- primary_split(n, if (k == 1) APP_SEED else APP_SEED + k); s <- sp$scr
    sc <- screen_gselo(co$data$time[s], co$data$status[s], co$X[s, , drop = FALSE], co$tau, cap = cap)
    sel[[k]] <- colnames(co$X)[sc$A]; sz[k] <- length(sc$A)
  }
  fr <- sort(table(unlist(sel)) / n_split, decreasing = TRUE)
  stem <- file.path(od, sprintf("%s_cap%d", cohort, cap))
  utils::write.csv(data.frame(var = names(fr), freq = as.numeric(fr)), paste0(stem, "_frequency.csv"), row.names = FALSE)
  sp <- primary_split(n, APP_SEED)
  p <- analyze_split(co$data$time, co$data$status, co$X, co$tau, sp, "gselo", cap = cap, subsample = FALSE)
  utils::write.csv(p$table, paste0(stem, "_primary.csv"), row.names = FALSE)
  at_cap <- mean(sz >= cap)
  message(sprintf("%s, cap %d: selected set at the cap in %.0f%% of %d splits", cohort, cap, 100 * at_cap, n_split))
  invisible(at_cap)
}
