## =====================================================================
## File 01_ah_core.R -- Lin-Ying additive hazards estimating equations with
## risk-set centering on [0, tau]:
##   b = n^-1 sum_i int Z_i(t) dN_i(t),  V = n^-1 sum_i int Y_i(t) Z_i Z_i^T dt,
##   Z_i(t) = X_i - Xbar(t).
## (a) ah_bV(): ambient-dimension b and V for screening, computed in row
##     blocks (memory O(B p + p^2)), exact in finite arithmetic.
## (b) ah_setup(): per-subject quantities on a low-dimensional set A for the
##     refit stage (full refit, IPW refits, sandwich components).
##
## Notation for (b): records sorted by observed time; groups k = 1..K of tied
## times t_k; S0_k, S1_k = risk-set count / covariate sum on [t_k, inf);
## Delta_k = length of (t_{k-1}, t_k] within [0, tau]; d_k = events at t_k
## (counted only if t_k <= tau).  Cumulatives up to group g:
##   Cc_g = sum_{l<=g} Delta_l Xbar_l,       Nb_g = sum_{l<=g} d_l / S0_l,
##   Eb_g = sum_{l<=g} Xbar_l d_l / S0_l,    Dc_g = sum_{l<=g} Delta_l Xbar_l Xbar_l^T.
## For subject i in group g (u_i = min(T_i, tau)):
##   e_i  = int Z_i dN_i                   = d_i (X_i - Xbar_g)
##   G_i  = int Y_i Z_i X_i^T dt           = (u_i X_i - Cc_g) X_i^T
##   L_i  = int Y_i Z_i Z_i^T dt           = u_i X_i X_i^T - X_i Cc_g^T - Cc_g X_i^T + Dc_g
##   q_i(beta)  = e_i - G_i beta                         (score contribution, eq. 4)
##   xi_i(beta) = e_i - (X_i Nb_g - Eb_g) - L_i beta     (Breslow martingale residual)
## =====================================================================

## Basic ordering / grouping on [0, tau]
ah_order <- function(time, status, tau) {
  o <- order(time)
  t <- time[o]; s <- status[o]
  first <- match(t, t)                        # first index of each tie group
  gstart <- which(!duplicated(t)); K <- length(gstart)
  grp <- cumsum(!duplicated(t))              # group id per sorted record
  tg <- t[gstart]
  tc <- pmin(tg, tau)
  Delta <- diff(c(0, tc))
  ev <- as.numeric(s == 1 & t <= tau)
  dk <- tabulate(grp[ev == 1], nbins = K)
  S0 <- length(t) - gstart + 1
  list(o = o, t = t, ev = ev, grp = grp, gstart = gstart, K = K, tg = tg,
       Delta = Delta, dk = dk, S0 = S0, u = pmin(t, tau), n = length(t), tau = tau)
}

## ---------------------------------------------------------------------
## (a) ambient b and V (screening), blockwise
## ---------------------------------------------------------------------
ah_bV <- function(time, status, X, tau, block = 5000L) {
  ob <- ah_order(time, status, tau)
  X <- X[ob$o, , drop = FALSE]
  n <- ob$n; p <- ncol(X)
  ## V = n^-1 [ sum_i u_i X_i X_i^T - sum_k Delta_k S1_k S1_k^T / S0_k ]
  ## b = n^-1 [ sum_i ev_i X_i - sum_k d_k S1_k / S0_k ]
  V <- crossprod(X * sqrt(ob$u))
  b <- colSums(X * ob$ev)
  carry <- numeric(p)                       # sum of rows after the current block
  starts <- rev(seq(1L, n, by = block))
  for (a in starts) {
    e <- min(n, a + block - 1L)
    Xb <- X[a:e, , drop = FALSE]
    rc <- apply(Xb, 2, function(x) rev(cumsum(rev(x))))
    if (is.null(dim(rc))) rc <- matrix(rc, nrow = 1)
    rc <- sweep(rc, 2, carry, "+")          # reverse cumsum at each row
    gs <- ob$gstart[ob$gstart >= a & ob$gstart <= e]
    if (length(gs)) {
      k <- ob$grp[gs]
      S1 <- rc[gs - a + 1L, , drop = FALSE]
      w1 <- ob$Delta[k] / ob$S0[k]
      keep <- w1 > 0
      if (any(keep)) V <- V - crossprod(S1[keep, , drop = FALSE] * sqrt(w1[keep]))
      w2 <- ob$dk[k] / ob$S0[k]
      if (any(w2 > 0)) b <- b - colSums(S1 * w2)
    }
    carry <- rc[1, ]
  }
  list(b = b / n, V = V / n, n = n, n_event = sum(ob$ev), diagV = diag(V) / n)
}

## ---------------------------------------------------------------------
## (b) low-dimensional set-up on A (columns of XA)
## ---------------------------------------------------------------------
ah_setup <- function(time, status, XA, tau, w = NULL) {
  ob <- ah_order(time, status, tau)
  XA <- as.matrix(XA)[ob$o, , drop = FALSE]
  d <- ncol(XA); K <- ob$K
  ww <- if (is.null(w)) rep(1, ob$n) else w[ob$o]       # analysis weights (survey sensitivity)
  rc <- apply(XA * ww, 2, function(x) rev(cumsum(rev(x))))
  if (is.null(dim(rc))) rc <- matrix(rc, ncol = d)
  S1 <- rc[ob$gstart, , drop = FALSE]
  S0 <- rev(cumsum(rev(ww)))[ob$gstart]
  dkw <- as.numeric(tapply(ww * ob$ev, factor(ob$grp, levels = seq_len(K)), sum)); dkw[is.na(dkw)] <- 0
  Xbar <- S1 / S0
  zero <- S0 <= 0                                   # no weighted risk (possible with zero weights)
  if (any(zero)) Xbar[zero, ] <- 0
  Cc <- apply(Xbar * ob$Delta, 2, cumsum); if (is.null(dim(Cc))) Cc <- matrix(Cc, ncol = d)
  dL <- ifelse(zero, 0, dkw / S0)
  Nb <- cumsum(dL)
  Eb <- apply(Xbar * dL, 2, cumsum); if (is.null(dim(Eb))) Eb <- matrix(Eb, ncol = d)
  g <- ob$grp
  e <- (XA - Xbar[g, , drop = FALSE]) * ob$ev
  P <- ob$u * XA - Cc[g, , drop = FALSE]               # G_i = P_i X_i^T
  R0 <- XA * Nb[g] - Eb[g, , drop = FALSE]             # baseline part of xi_i
  W <- sum(ww)
  ## full V = W^-1 [sum w u X X^T - sum_k Delta_k S0_k Xbar Xbar^T]
  Vfull <- (crossprod(XA * sqrt(ww * ob$u)) - crossprod(Xbar * sqrt(ob$Delta * S0))) / W
  structure(list(ob = ob, X = XA, d = d, n = ob$n, g = g, Xbar = Xbar, Cc = Cc,
                 Nb = Nb, Eb = Eb, e = e, P = P, R0 = R0, Vfull = Vfull, w = ww, W = W,
                 bfull = colSums(e * ww) / W, ev = ob$ev, u = ob$u),
            class = "ah_setup")
}

## Dc_g beta for every subject (n x d) given beta
.Dc_beta <- function(st, beta) {
  v <- st$Xbar * (st$ob$Delta * drop(st$Xbar %*% beta))
  cs <- apply(v, 2, cumsum); if (is.null(dim(cs))) cs <- matrix(cs, ncol = st$d)
  cs[st$g, , drop = FALSE]
}
## per-subject q_i(beta) (n x d) and xi_i(beta) (n x d), in sorted order
q_scores <- function(st, beta) st$e - st$P * drop(st$X %*% beta)
xi_scores <- function(st, beta) {
  xb <- drop(st$X %*% beta)
  Lb <- st$u * st$X * xb - st$X * drop(st$Cc[st$g, , drop = FALSE] %*% beta) -
        st$Cc[st$g, , drop = FALSE] * xb + .Dc_beta(st, beta)
  st$e - st$R0 - Lb
}

## Full (split-full) refit with martingale sandwich
ah_fit_full <- function(st, ridge = 1e-8) {
  beta <- solve_ridge(st$Vfull, st$bfull, ridge)
  xi <- xi_scores(st, beta)
  Om <- crossprod(xi) / st$n
  Hi <- solve(st$Vfull + diag(ridge, st$d))
  Sig <- Hi %*% (Om / st$n) %*% t(Hi)
  list(beta = beta, Sigma = Sig, se = sqrt(diag(Sig)), Omega = Om, H = st$Vfull)
}

## Sum over sampled records (with multiplicities cnt, sorted order) of
## weights * L_i, computed without storing Dc for every group.
.sum_wL <- function(st, w) {
  idx <- which(w != 0); ww <- w[idx]
  Xi <- st$X[idx, , drop = FALSE]; Ci <- st$Cc[st$g[idx], , drop = FALSE]
  A1 <- crossprod(Xi * (ww * st$u[idx]), Xi) - crossprod(Xi * ww, Ci) - crossprod(Ci * ww, Xi)
  ## sum_i w_i Dc_{g_i} = sum_l Delta_l Xbar_l Xbar_l^T * W_l, W_l = sum_{i: g_i >= l} w_i
  Wg <- as.numeric(tapply(ww, factor(st$g[idx], levels = seq_len(st$ob$K)), sum))
  Wg[is.na(Wg)] <- 0
  Wl <- rev(cumsum(rev(Wg)))
  A1 + crossprod(st$Xbar * sqrt(pmax(st$ob$Delta * Wl, 0)))
}
