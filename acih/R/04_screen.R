## =====================================================================
## File 04_screen.R -- screening on D_scr (Section 2.1, Web Appendix B):
##  * GSELO continuation path for
##      Q(beta) = beta^T V beta / 2 - b^T beta + lambda* sum_j {1 - exp(-theta beta_j^2)},
##    lambda* = log{max(N_event, 1)} / n_scr, ascending theta grid with warm
##    starts, MM (local quadratic) updates
##      beta <- (V + diag(2 lambda* theta exp(-theta beta^2)) + ridge I)^-1 b,
##    at most 800 iterations, un-ridged KKT tolerance 1e-3 (sup norm);
##  * size-constrained joint EBIC search over the theta x threshold grid,
##    with ridge-stabilised unpenalised refits on the candidate supports and
##    criterion 2 n_scr * loss / phi + |A| log D_scr + 2 gamma log choose(p, |A|)
##    (D_scr = number of screening-split events; log n_scr with ebic_n = "n");
##  * marginal pre-screen u_j = |b_j| / sqrt(V_jj), top K;
##  * fixed-budget ranking (real data): top-K absolute coefficients of the
##    unpenalised multivariable fit.
## =====================================================================

theta_grid_default <- function()
  c(0.1, 0.25, 0.5, seq(1, 25, by = 2), seq(30, 150, by = 5), seq(160, 300, by = 20),
    seq(350, 500, by = 50), seq(600, 1000, by = 100), 1200, 1500, 2000, 2500, 3000, 4000, 5000)
eta_grid_default <- function()
  c(0, 0.0002, 0.0005, 0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.015, 0.02,
    0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, Inf)

gselo_kkt <- function(V, b, beta, lam, theta)
  max(abs(drop(V %*% beta) - b + 2 * lam * theta * beta * exp(-theta * beta^2)))

gselo_path <- function(V, b, lam, thetas = theta_grid_default(), maxit = 800,
                       tol = 1e-3, ridge = 1e-8) {
  p <- length(b); beta <- numeric(p)
  path <- matrix(0, p, length(thetas)); kkt <- numeric(length(thetas)); it_used <- integer(length(thetas))
  conv <- logical(length(thetas))
  for (k in seq_along(thetas)) {
    th <- thetas[k]
    for (it in seq_len(maxit)) {
      D <- 2 * lam * th * exp(-th * beta^2)
      A <- V; diag(A) <- diag(A) + D + ridge
      R <- tryCatch(chol(A), error = function(e) NULL)
      bn <- if (is.null(R)) solve_ridge(A, b, 1e-8) else backsolve(R, forwardsolve(t(R), b))
      dlt <- max(abs(bn - beta)); beta <- bn
      kk <- gselo_kkt(V, b, beta, lam, th)
      if (kk <= tol) break
      if (dlt < 1e-12) break
    }
    path[, k] <- beta; kkt[k] <- kk; it_used[k] <- it; conv[k] <- kk <= tol
  }
  list(path = path, theta = thetas, kkt = kkt, iter = it_used, converged = conv, lambda = lam)
}

## ridge-stabilised unpenalised refit on A and its un-ridged loss
refit_loss <- function(V, b, A, ridge = 1e-8) {
  if (!length(A)) return(0)
  bA <- solve_ridge(V[A, A, drop = FALSE], b[A], ridge)
  0.5 * sum(bA * (V[A, A, drop = FALSE] %*% bA)) - sum(b[A] * bA)
}

ebic_select <- function(V, b, gp, n_scr, gamma = 0.25, etas = eta_grid_default(),
                        cap = 20, p_universe = length(b), ridge = 1e-8, loss_mult = 1,
                        log_pen = log(n_scr)) {
  cache <- new.env(); best <- list(crit = Inf)
  tab <- vector("list", length(gp$theta) * length(etas)); kk <- 0
  for (k in seq_along(gp$theta)) for (eta in etas) {
    A <- which(abs(gp$path[, k]) > eta)
    if (length(A) > cap) next
    key <- paste0("A:", paste(A, collapse = ","))
    loss <- cache[[key]]
    if (is.null(loss)) { loss <- refit_loss(V, b, A, ridge); assign(key, loss, envir = cache) }
    crit <- loss_mult * n_scr * loss + length(A) * log_pen + 2 * gamma * lchoose(p_universe, length(A))
    kk <- kk + 1; tab[[kk]] <- c(theta = gp$theta[k], eta = eta, size = length(A), crit = crit)
    if (crit < best$crit - 1e-12 || (abs(crit - best$crit) <= 1e-12 && length(A) < length(best$A)))
      best <- list(crit = crit, A = A, theta = gp$theta[k], eta = eta, k = k)
  }
  best$table <- do.call(rbind, tab[seq_len(kk)])
  ## the criterion depends on (theta, eta) only through the selected set, so
  ## the same set is usually produced on a range of theta values; a boundary
  ## selection means that the chosen set is produced only at an endpoint of
  ## the theta grid
  th_best <- unique(unlist(lapply(seq_along(gp$theta), function(k)
    if (any(vapply(etas, function(eta) setequal(which(abs(gp$path[, k]) > eta), best$A), TRUE))) gp$theta[k])))
  best$theta_range <- range(th_best)
  best$boundary_low <- all(th_best == min(gp$theta))
  best$boundary_high <- all(th_best == max(gp$theta))
  best$boundary <- best$boundary_low || best$boundary_high
  best$kkt <- gp$kkt[best$k]; best$converged <- gp$converged[best$k]
  best
}

## Complete GSELO screen on D_scr
## time_scale = "rate" (default): the screen is computed on the time scale on
## which the crude event rate D / T on [0, tau] of the screening data equals
## one, and the EBIC uses the quasi-likelihood scaling 2 n loss.  The
## additive-hazards loss l(beta) = beta' V beta / 2 - b' beta scales with the
## time unit while the EBIC and GSELO penalties do not; on the unit-rate scale
## 2 n {l(A without j) - l(A)} approximates the Wald chi-square of coordinate j
## (the predictable variation of the score is close to V when the hazard is
## close to its average), so the criterion is
##   2 n l_A / phi + |A| log n + 2 gamma log choose(p, |A|),  phi = D / T,
## and the selected set does not depend on the unit of time.  The refit and
## the inference are equivariant under rescaling of time and are unaffected.
## time_scale = "none" reproduces the unscaled criterion n l_A + ... of
## acih <= 1.0.2.
time_rate <- function(time, status, tau) sum(status == 1 & time <= tau) / sum(pmin(time, tau))

screen_gselo <- function(time, status, X, tau, gamma = 0.25, thetas = theta_grid_default(),
                         etas = eta_grid_default(), cap = 20, candidates = NULL,
                         time_scale = c("rate", "none"), ebic_n = getOption("acih.ebic_n", "events")) {
  ## ebic_n: effective sample size in the BIC term |A| log(.) of the EBIC:
  ## "events" (default from acih 1.3.1; number of events on [0, tau], as for
  ## survival BIC, Volinsky and Raftery 2000) or "n" (screening-split size,
  ## the criterion of acih <= 1.3.0)
  t0 <- tic()
  time_scale <- match.arg(time_scale); ebic_n <- match.arg(ebic_n, c("n", "events"))
  sc_fac <- if (time_scale == "rate") time_rate(time, status, tau) else 1
  time <- time * sc_fac; tau <- tau * sc_fac
  Xs <- if (is.null(candidates)) X else X[, candidates, drop = FALSE]
  bv <- ah_bV(time, status, Xs, tau)
  lam <- log(max(bv$n_event, 1)) / bv$n
  gp <- gselo_path(bv$V, bv$b, lam, thetas)
  sel <- ebic_select(bv$V, bv$b, gp, bv$n, gamma, etas, cap,
                     p_universe = if (is.null(candidates)) ncol(X) else length(candidates),
                     loss_mult = if (time_scale == "rate") 2 else 1,
                     log_pen = if (ebic_n == "events") log(max(bv$n_event, 2)) else log(bv$n))
  A <- if (is.null(candidates)) sel$A else candidates[sel$A]
  list(A = A, theta = sel$theta, theta_range = sel$theta_range, eta = sel$eta, boundary_low = sel$boundary_low,
       boundary_high = sel$boundary_high, kkt = sel$kkt, converged = sel$converged,
       all_converged = all(gp$converged), lambda = lam, time = tic() - t0,
       min_eig_V = NA_real_, time_factor = sc_fac, ebic_n = ebic_n)
}

## Marginal pre-screen (direct top-K) with optional candidate-restricted GSELO
screen_marginal <- function(time, status, X, tau, K, then_gselo = FALSE, ...) {
  t0 <- tic()
  bv <- ah_bV(time, status, X, tau)
  u <- abs(bv$b) / sqrt(pmax(bv$diagV, 1e-300))
  cand <- order(-u, seq_along(u))[seq_len(K)]
  out <- list(A = sort(cand), candidates = sort(cand), time = tic() - t0)
  if (then_gselo) {
    g <- screen_gselo(time, status, X, tau, candidates = sort(cand), ...)
    out$A <- g$A; out$gselo <- g; out$time <- tic() - t0
  }
  out
}

## Fixed-budget ranking used in the applications: top-K |coefficients| of
## the unpenalised multivariable fit on D_scr
screen_topk_unpen <- function(time, status, X, tau, K) {
  bv <- ah_bV(time, status, X, tau)
  beta <- solve_ridge(bv$V, bv$b, 1e-8)
  list(A = sort(order(-abs(beta))[seq_len(K)]), beta = beta)
}
