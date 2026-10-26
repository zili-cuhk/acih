## =====================================================================
## File 07_extensions.R -- Poisson (Bernoulli) subsampling design.
##
## Each record of D_inf is included independently with probability p_i,
## fixed before sampling, with sum_i p_i = r (expected subsample size).
##   b^P = m^-1 sum_i R_i / p_i e_i,     V^P = m^-1 sum_i R_i / p_i G_i,
##   Var*(score) = m^-2 sum_i (1 - p_i) / p_i q_i q_i^T            (exact),
##   estimated by the Horvitz-Thompson form m^-2 sum_i R_i (1 - p_i) / p_i^2 q_i q_i^T.
## Optimal design for tr(C H^-1 Var* H^-T C^T): p_i = min(1, kappa sqrt(a_i)).
## Event-complete design: p_i = 1 for observed events (when #events < r), the
## optimal (or uniform) design among censored records.  Stabilisation mixes
## the normalised scores with the uniform design: s~_i = (1-eps) s_i/sum s + eps/m.
## =====================================================================

## capped proportional allocation with total r
.capped_alloc <- function(w, total) {
  p <- numeric(length(w)); free <- rep(TRUE, length(w))
  for (it in 1:200) {
    budget <- total - sum(!free)
    if (budget <= 0) break
    p[free] <- budget * w[free] / sum(w[free])
    over <- free & p >= 1
    if (!any(over)) break
    p[over] <- 1; free <- free & !over
  }
  pmin(p, 1)
}

## s: nonnegative scores (n), ev: event indicators, r: expected size
poisson_inclusion <- function(s, ev, r, eps = 0.10, event_complete = TRUE) {
  n <- length(s)
  if (r >= n) return(list(p = rep(1, n), event_complete = TRUE))
  stab <- function(x) { x <- if (sum(x) > 0) x / sum(x) else rep(1 / length(x), length(x)); (1 - eps) * x + eps / length(x) }
  ne <- sum(ev)
  if (event_complete && ne < r) {
    p <- numeric(n); p[ev == 1] <- 1
    cen <- which(ev == 0)
    p[cen] <- .capped_alloc(stab(s[cen]), r - ne)
    return(list(p = p, event_complete = TRUE))
  }
  list(p = .capped_alloc(stab(s), r), event_complete = FALSE)
}
poisson_probs <- function(s, ev, r, eps = 0.10) poisson_inclusion(s, ev, r, eps, TRUE)$p

## allocation scores on the inference split from the screening-split pilot
poisson_scores <- function(st, pilot, rule = "inference", C = NULL) {
  if (rule == "uniform") return(rep(1, st$n))
  q <- q_scores(st, pilot$beta)
  if (is.null(C)) C <- diag(st$d)
  sqrt(rowSums((q %*% t(C %*% pilot$Hinv))^2))
}

## Poisson IPW refit with total sandwich; Up: uniforms (length n, sorted order)
ah_fit_poisson <- function(st, p, Up, ridge = 1e-8) {
  inc <- Up < p
  w <- inc / p / st$n
  b <- colSums(st$e * w); H <- crossprod(st$P * w, st$X)
  beta <- solve_ridge(H, b, ridge)
  xi <- xi_scores(st, beta); Om <- crossprod(xi) / st$n
  q <- q_scores(st, beta)
  Gp <- crossprod(q * sqrt(inc * (1 - p) / p^2)) / st$n^2     # HT estimate of Var*(score)
  Hi <- solve(H + diag(ridge, st$d))
  Sig <- Hi %*% (Om / st$n + Gp) %*% t(Hi)
  list(beta = beta, se = sqrt(pmax(diag(Sig), 0)), Sigma = Sig,
       Sigma_omega_only = Hi %*% (Om / st$n) %*% t(Hi),
       Sigma_gamma_only = Hi %*% Gp %*% t(Hi),
       gamma_frac = mean(diag(Hi %*% Gp %*% t(Hi)) / pmax(diag(Sig), 1e-300)),
       size = sum(inc), n_event_sampled = sum(inc * st$ev), max_w = max(1 / p))
}

## ---------------------------------------------------------------------
## Case-cohort comparator: Kulich and Lin (2000, Biometrika 87, 73-87),
## additive hazards regression for case-cohort studies, Bernoulli design.
## All cases (events in [0, tau]) are included with weight 1 and controls
## with probability p_i, weight 1/p_i:  r_i = d_i + (1 - d_i) j_i / p_i.
## Risk-set averages use the sampled records only (covariates are not
## available outside the case-cohort sample):
##   Zbar_H(t) = sum r_j Y_j Z_j / sum r_j Y_j,
##   b_H = [sum r_i int (Z_i - Zbar_H)^2 Y_i dt]^-1 sum int (Z_i - Zbar_H) dN_i,
## variance n^-1 D^-1 (S_A + S_H) D^-1 with
##   S_A = n^-1 sum int (Z_i - Zbar_H)^2 dN_i,  D = n^-1 sum r_i int (Z_i - Zbar_H)^2 Y_i dt,
##   S_H = n^-1 sum (1 - p_i) / p_i^2 j_i (1 - d_i) S_i^2,
##   S_i = int (Z_i - Zbar_H) {dN_i - Y_i dL_H - b_H' Z_i Y_i dt}  (their eqs. 3, 6, 7).
## st: ah_setup of the inference split (sorted order); p: control inclusion
## probabilities and Up: uniforms, both in sorted order (shared with PUIPW).
## ---------------------------------------------------------------------
ah_fit_casecohort <- function(st, p, Up, ridge = 1e-8) {
  inc <- (Up < p) | st$ev == 1
  r <- ifelse(st$ev == 1, 1, (Up < p) / p)
  sc <- ah_setup(st$ob$t, st$ev, st$X, st$ob$tau, w = r)     # input already sorted
  n <- st$n; d <- st$d
  beta <- solve_ridge(sc$Vfull, sc$bfull, ridge)
  D <- sc$Vfull * sc$W / n
  SA <- crossprod(sc$e) / n                                   # e_i = 0 for non-cases
  Si <- xi_scores(sc, beta)
  wH <- ifelse(sc$ev == 0 & inc, (1 - p) / p^2, 0)
  SH <- crossprod(Si * sqrt(wH)) / n
  Di <- solve(D + diag(ridge, d))
  Sig <- Di %*% (SA + SH) %*% t(Di) / n
  VH <- Di %*% SH %*% t(Di) / n
  list(beta = beta, se = sqrt(pmax(diag(Sig), 0)), Sigma = Sig,
       gamma_frac = mean(diag(VH) / pmax(diag(Sig), 1e-300)),
       size = sum(inc), max_w = max(r))
}
