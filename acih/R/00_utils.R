## =====================================================================
## acih: Active-coordinate inference for sparse additive hazards regression
## under informative subsampling.   File 00_utils.R -- helpers (base R).
## =====================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
tic <- function() proc.time()[["elapsed"]]

ar1_cor <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))

## AR(1) Gaussian rows without forming a p x p Cholesky (O(np))
rar1 <- function(n, p, rho) {
  Z <- matrix(stats::rnorm(n * p), n, p)
  if (rho == 0) return(Z)
  s <- sqrt(1 - rho^2)
  for (j in 2:p) Z[, j] <- rho * Z[, j - 1] + s * Z[, j]
  Z
}
## Block-exchangeable Gaussian rows: blocks of size `bs` with within-block
## correlation rho, independent across blocks
rblock <- function(n, p, rho, bs = 10) {
  nb <- ceiling(p / bs)
  F <- matrix(stats::rnorm(n * nb), n, nb)[, rep(seq_len(nb), each = bs)[seq_len(p)], drop = FALSE]
  sqrt(rho) * F + sqrt(1 - rho) * matrix(stats::rnorm(n * p), n, p)
}

## Solve with a fixed small ridge (the "numerical safeguard" of Section 2.3)
solve_ridge <- function(A, b, ridge = 1e-8) {
  A <- A + diag(ridge, nrow(A))
  out <- tryCatch(solve(A, b), error = function(e) NULL)
  if (is.null(out)) out <- qr.solve(A + diag(1e-6 * max(1, mean(abs(diag(A)))), nrow(A)), b)
  out
}

wilson <- function(x, n, z = 1.96) {
  p <- x / n; den <- 1 + z^2 / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den)
}

msg <- function(...) message(sprintf(...))
