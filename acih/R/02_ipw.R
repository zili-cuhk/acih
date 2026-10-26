## =====================================================================
## File 02_ipw.R -- allocation rules, with-replacement IPW refits and the
## total sandwich covariance (Sections 2.2-2.4).
##
## Score types:
##   "q"  : the manuscript's IPW equation b*_A - V*_A beta, built from
##          q_i(beta) = int Z_i {dN_i - Y_i X_i^T beta dt}           (eq. 7)
##   "xi" : residual (control-variate) IPW equation built from the Breslow
##          martingale residual xi_i(beta); it has the same full-split total
##          as "q" (sum_i Y_i Z_i(t) = 0), so it is conditionally unbiased for
##          the same full-inference equation, its derivative is symmetric,
##          and its subsampling variance omits the baseline-hazard term that
##          dominates Gamma for "q".  (Suggested extension; see review.)
## =====================================================================

## ---- allocation ---------------------------------------------------------
## scores: n x d matrix of pilot scores (sorted order); Hinv: pilot inverse
## information; C: contrast matrix (rows); type:
##   "uniform", "inference" (sqrt of a_i = ||C H^-1 s_i||^2, eq. 11),
##   "scorenorm" (||s_i||), "zuo" (see alloc_zuo()).
## Stratified by event status with pre-stabilisation masses omega, then
## stabilised: pi = (1 - eps) pi_str + eps / n.
allocate <- function(scores, ev, type = "inference", Hinv = NULL, C = NULL,
                     eps = 0.10, omega = c(0.5, 0.5), zuo_extra = NULL) {
  n <- nrow(scores)
  if (type == "uniform") return(rep(1 / n, n))
  sz <- switch(type,
    inference = {
      if (is.null(C)) C <- diag(ncol(scores))
      M <- scores %*% t(C %*% Hinv)          # rows: C H^-1 s_i
      sqrt(rowSums(M^2))
    },
    scorenorm = sqrt(rowSums(scores^2)),
    zuo = alloc_zuo(scores, zuo_extra),
    stop("unknown allocation"))
  strata <- list(which(ev == 0), which(ev == 1))
  nonempty <- sapply(strata, length) > 0
  om <- omega * nonempty; om <- om / sum(om)          # redistribute empty-stratum mass
  pis <- numeric(n)
  for (h in 1:2) {
    idx <- strata[[h]]; if (!length(idx)) next
    s <- sz[idx]
    pis[idx] <- if (sum(s) > 0) om[h] * s / sum(s) else om[h] / length(idx)
  }
  (1 - eps) * pis + eps / n
}

## Allocation of Zuo et al. (2021, Statist. Med.), eq. (11), "upper bound"
## form.  !!! VERIFY against eq. (11) of the paper before use: implemented
## here as the triangle-inequality bound ||e_i|| + ||G_i beta_pil|| of the
## score norm, using full-split centering, with the same event
## stratification and stabilisation as the other rules.
alloc_zuo <- function(scores, extra) {
  if (is.null(extra)) return(sqrt(rowSums(scores^2)))
  sqrt(rowSums(extra$e^2)) + sqrt(rowSums(extra$Pb^2))
}

## with-replacement draws by inverse CDF from common uniforms U (length r)
draw_ipw <- function(pi, U) {
  cp <- cumsum(pi); cp[length(cp)] <- 1
  pmin(findInterval(U, cp) + 1L, length(pi))
}

## ---- IPW refit and total sandwich ---------------------------------------
## st: ah_setup on the inference split; pi: probabilities (sorted order);
## I: sampled indices (sorted order); score: "q" or "xi".
ah_fit_ipw <- function(st, pi, I, score = c("q", "xi"), ridge = 1e-8,
                       omega_only = FALSE) {
  score <- match.arg(score)
  n <- st$n; r <- length(I); d <- st$d
  cnt <- tabulate(I, nbins = n)
  w <- cnt / (n * pi) / r                 # (1/r) sum_k 1/(n pi_{I_k}) per record
  if (score == "q") {
    b <- colSums(st$e * w)
    H <- crossprod(st$P * w, st$X)        # sum w_i (u X - Cc) X^T  (not symmetric)
  } else {
    b <- colSums((st$e - st$R0) * w)
    H <- .sum_wL(st, w)
  }
  beta <- solve_ridge(H, b, ridge)
  xi <- xi_scores(st, beta)
  Om <- crossprod(xi) / n
  s <- if (score == "q") q_scores(st, beta) else xi
  sbar <- colMeans(s)
  Ga <- crossprod(s / sqrt(pi)) / n^2 - tcrossprod(sbar)
  Hi <- solve(H + diag(ridge, d))
  M <- Om / n + (if (omega_only) 0 else Ga / r)
  Sig <- Hi %*% M %*% t(Hi)
  list(beta = beta, Sigma = Sig, se = sqrt(pmax(diag(Sig), 0)), Omega = Om, Gamma = Ga,
       H = H, gamma_frac = mean(diag(Hi %*% (Ga / r) %*% t(Hi)) / pmax(diag(Sig), 1e-300)),
       Sigma_omega_only = Hi %*% (Om / n) %*% t(Hi),
       Sigma_gamma_only = Hi %*% (Ga / r) %*% t(Hi),
       ess = sum(cnt)^2 / sum(cnt^2), n_event_sampled = sum(cnt * st$ev),
       max_w = max(1 / (n * pi)), min_pi_n = min(pi) * n)
}

## Pilot quantities on the screening split for allocation
pilot_fit <- function(time, status, XA, tau, ridge = 1e-8) {
  st <- ah_setup(time, status, XA, tau)
  beta <- solve_ridge(st$Vfull, st$bfull, ridge)
  list(beta = beta, H = st$Vfull, Hinv = solve(st$Vfull + diag(ridge, st$d)))
}

## Probabilities for a rule on the inference split (sorted order)
ipw_probs <- function(st, pilot, rule, score = "q", C = NULL, eps = 0.10,
                      omega = c(0.5, 0.5)) {
  if (rule == "uniform") return(rep(1 / st$n, st$n))
  s <- if (score == "q") q_scores(st, pilot$beta) else xi_scores(st, pilot$beta)
  extra <- if (rule == "zuo") list(e = st$e, Pb = st$P * drop(st$X %*% pilot$beta)) else NULL
  allocate(s, st$ev, rule, Hinv = pilot$Hinv, C = C, eps = eps, omega = omega, zuo_extra = extra)
}
