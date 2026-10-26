## =====================================================================
## File 03_datagen.R -- data generation under the additive hazards model
##   lambda(t | X) = lambda0(t) + X^T beta0,
## covariates bounded by `bound` (componentwise clipping), baseline chosen
## before simulation so that the hazard is nonnegative for every X in the
## support, censoring calibrated once per design.
## =====================================================================

## Covariates: "ar1" Gaussian, "block" exchangeable blocks, "t5" (AR(1)
## Gaussian copula with t5 margins scaled to unit variance), each clipped at
## +-bound.  `binary_every` > 0 dichotomises every k-th covariate (0/1).
gen_X <- function(n, p, design = "ar1", rho = 0.5, bound = 2.5, binary_every = 0,
                  block_size = 10) {
  X <- switch(design,
    ar1 = rar1(n, p, rho),
    block = rblock(n, p, rho, block_size),
    t5 = { Z <- rar1(n, p, rho); stats::qt(stats::pnorm(Z), 5) / sqrt(5 / 3) },
    stop("unknown covariate design"))
  X <- pmin(pmax(X, -bound), bound)
  if (binary_every > 0) {
    jb <- which(seq_len(p) %% binary_every == 0)
    X[, jb] <- (X[, jb] > 0) * 1
  }
  X
}

## Baseline: list(type = "constant", level) or list(type = "linear",
## level, slope): lambda0(t) = level + slope * t.  The level is set from
## the covariate bound: level = bound * ||beta0||_1 + margin.
make_baseline <- function(beta0, bound = 2.5, margin = 0.5, type = "constant", slope_frac = 0) {
  lev <- bound * sum(abs(beta0)) + margin
  list(type = type, level = lev, slope = slope_frac * lev)
}

## Event times by inversion: Lambda(t|X) = level t + slope t^2 / 2 + (X^T b) t
gen_T <- function(X, beta0, base, tv = NULL) {
  n <- nrow(X); E <- stats::rexp(n)
  a <- base$level + drop(X %*% beta0)
  if (any(a < 0)) a <- pmax(a, 0)          # only reachable under misspecification designs
  if (!is.null(tv)) {
    ## time-varying effect of covariate tv$j: beta_j(t) = beta_j + tv$slope * (t - tv$center)
    a <- a + tv$slope * (-tv$center) * X[, tv$j]
    c2 <- base$slope / 2 + tv$slope * X[, tv$j] / 2
  } else c2 <- rep(base$slope / 2, n)
  ## solve c2 t^2 + a t = E  (c2 may be 0)
  t <- ifelse(abs(c2) < 1e-12, E / a, (-a + sqrt(pmax(a^2 + 4 * c2 * E, 0))) / (2 * c2))
  t[!is.finite(t) | t < 0] <- Inf
  t
}

## Censoring: "uniform" C ~ U(0, cmax); "covexp" C ~ Exp(c0 exp(g X1));
## "weibull" shape 0.8, scale c0; "frailty" exponential censoring sharing a
## gamma frailty with the event time (violates conditional independence).
gen_C <- function(X, cens, par, frail = NULL) {
  n <- nrow(X)
  switch(cens,
    uniform = stats::runif(n, 0, par),
    covexp = stats::rexp(n, par * exp(0.6 * X[, 1])),
    weibull = stats::rweibull(n, 0.8, par),
    frailty = stats::rexp(n, par * frail),
    none = rep(Inf, n),
    stop("unknown censoring"))
}

## Scenario object
make_design <- function(n = 100000, p = 500, s_design = "s4", rho = 0.5,
                        xdesign = "ar1", cens = "uniform", cens_rate = 0.5,
                        tau = 0.3, base_type = "constant", slope_frac = 0,
                        bound = 2.5, binary_every = 0, tv = NULL,
                        r = c(1000, 5000), label = NULL) {
  beta0 <- numeric(p)
  if (s_design == "s4") beta0[1:4] <- c(0.35, -0.35, 0.32, -0.32)
  else if (s_design == "s10") beta0[1:10] <- c(0.50, -0.50, 0.40, -0.40, 0.30, -0.30, 0.25, -0.25, 0.20, -0.20)
  else if (s_design == "weak4") beta0[1:4] <- c(0.20, -0.20, 0.18, -0.18)
  else if (s_design == "cancel4") beta0[1:4] <- c(0.35, 0.35, 0.32, 0.32)   # used with negative within-pair correlation
  else if (s_design == "reinf4") beta0[1:4] <- c(0.70, 0.70, 0.65, 0.65)
  else if (s_design == "alt4") beta0[1:4] <- c(0.70, -0.70, 0.65, -0.65)
  else stop("unknown s_design")
  list(n = n, p = p, beta0 = beta0, S0 = which(beta0 != 0), rho = rho, xdesign = xdesign,
       cens = cens, cens_rate = cens_rate, cens_par = NA_real_, tau = tau,
       base = make_baseline(beta0, bound, type = base_type, slope_frac = slope_frac),
       bound = bound, binary_every = binary_every, tv = tv, r = r,
       label = label %||% sprintf("n%d_p%d_%s_%s_%s%.0f", n, p, s_design, xdesign, cens, 100 * cens_rate))
}

.gen_raw <- function(des, n, frail_var = 0) {
  X <- gen_X(n, des$p, des$xdesign, des$rho, des$bound, des$binary_every)
  frail <- if (des$cens == "frailty") stats::rgamma(n, 1 / 0.5, 1 / 0.5) else NULL
  if (!is.null(frail)) {
    ## shared frailty acts additively on the event hazard: extra rate (frail - 1) * base level / 2
    T <- gen_T(X, des$beta0, des$base, des$tv)
    rate2 <- pmax(frail - 1, 0) * des$base$level / 2    # zero rate means no extra events
    T2 <- rep(Inf, n); pos <- rate2 > 0
    T2[pos] <- stats::rexp(sum(pos), rate2[pos])
    T <- pmin(T, T2)
  } else T <- gen_T(X, des$beta0, des$base, des$tv)
  list(X = X, T = T, frail = frail)
}

## Calibrate the censoring parameter to the target proportion of censored
## observations (delta = 0) on an independent Monte Carlo population.
calibrate_design <- function(des, n_cal = 2e5, seed = 99) {
  set.seed(seed)
  if (des$cens == "none") {
    d2 <- des; d2$p <- max(des$S0, 2); d2$beta0 <- des$beta0[seq_len(d2$p)]
    raw <- .gen_raw(d2, n_cal)
    des$cens_par <- Inf; des$achieved_cens <- 0
    des$risk_frac_tau <- mean(raw$T >= des$tau); des$event_frac_tau <- mean(raw$T <= des$tau)
    return(des)
  }
  ## only the first max(S0) covariates enter the event and censoring models
  d2 <- des; d2$p <- max(des$S0, 2, if (!is.null(des$tv)) des$tv$j); d2$beta0 <- des$beta0[seq_len(d2$p)]
  raw <- .gen_raw(d2, n_cal)
  f <- function(lp) {
    C <- gen_C(raw$X[, 1:min(2, des$p), drop = FALSE], des$cens, exp(lp), raw$frail)
    mean(C < raw$T) - des$cens_rate
  }
  rng <- if (des$cens == "uniform") c(log(1e-4), log(1e3)) else c(log(1e-6), log(1e4))
  des$cens_par <- exp(stats::uniroot(f, rng, tol = 1e-10)$root)
  C <- gen_C(raw$X[, 1:min(2, des$p), drop = FALSE], des$cens, des$cens_par, raw$frail)
  time <- pmin(raw$T, C)
  des$achieved_cens <- mean(C < raw$T)
  des$risk_frac_tau <- mean(time >= des$tau)
  des$event_frac_tau <- mean(raw$T <= C & raw$T <= des$tau)
  des
}

gen_data <- function(des, n = des$n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  raw <- .gen_raw(des, n)
  C <- gen_C(raw$X[, 1:min(2, des$p), drop = FALSE], des$cens, des$cens_par, raw$frail)
  list(time = pmin(raw$T, C), status = as.integer(raw$T <= C), X = raw$X)
}

## Least-false (population) target on a fixed set A from a large sample:
## beta*_A = V_A^-1 b_A of the full estimating equation.
population_target <- function(des, A, n_big = 1e6, seed = 7) {
  ## generate only the first max(A) covariates (AR(1)/block margins do not
  ## depend on p); the event model uses the first p_gen coefficients
  set.seed(seed)
  pg <- max(A, des$S0); d2 <- des; d2$p <- pg; d2$beta0 <- des$beta0[seq_len(pg)]
  d <- gen_data(d2, n = n_big)
  st <- ah_setup(d$time, d$status, d$X[, A, drop = FALSE], des$tau)
  drop(solve(st$Vfull, st$bfull))
}
