## Package tests (run by R CMD check; a few minutes)
## 1. exactness of the blockwise ambient b, V; 2. score identities;
## 3. conditional unbiasedness of the IPW pieces; 4. one replicate of each
##    study type at reduced size.
library(acih)
.sum_wL <- acih:::.sum_wL; tic <- acih:::tic
set.seed(1); ok <- TRUE
chk <- function(name, val, tol) { pass <- is.finite(val) && val < tol; ok <<- ok && pass
  cat(sprintf("%-45s %-10.3g %s\n", name, val, if (pass) "PASS" else "FAIL")) }
des <- calibrate_design(make_design(n = 600, p = 6, cens = "uniform", cens_rate = 0.5), n_cal = 2e4)
d <- gen_data(des, seed = 2); d$time <- round(d$time, 2)
bv <- ah_bV(d$time, d$status, d$X, des$tau, block = 37); st <- ah_setup(d$time, d$status, d$X, des$tau)
chk("blockwise V equals in-memory V", max(abs(bv$V - st$Vfull)), 1e-12)
chk("blockwise b equals in-memory b", max(abs(bv$b - st$bfull)), 1e-12)
beta <- c(0.3, -0.2, 0.1, 0, 0.05, -0.1)
chk("mean q_i = b - V beta", max(abs(colMeans(q_scores(st, beta)) - (st$bfull - st$Vfull %*% beta))), 1e-12)
chk("mean xi_i = b - V beta", max(abs(colMeans(xi_scores(st, beta)) - (st$bfull - st$Vfull %*% beta))), 1e-12)
chk("sum_i L_i = n V", max(abs(.sum_wL(st, rep(1, st$n)) - st$n * st$Vfull)), 1e-9)
pil <- list(beta = beta, Hinv = solve(st$Vfull)); pi <- ipw_probs(st, pil, "inference")
acc <- 0; B <- 20000
for (k in 1:B) { I <- draw_ipw(pi, runif(30)); w <- tabulate(I, st$n) / (st$n * pi) / 30; acc <- acc + crossprod(st$P * w, st$X) }
chk("E*(V*) = V_full (MC, relative)", max(abs(acc / B - st$Vfull)) / max(abs(st$Vfull)), 0.02)
## exact conditional variance: Var*(score) = Gamma / r
q <- q_scores(st, beta); G <- crossprod(q / sqrt(pi)) / st$n^2 - tcrossprod(colMeans(q))
sims <- t(replicate(20000, { I <- draw_ipw(pi, runif(30)); colSums(q * (tabulate(I, st$n) / (st$n * pi) / 30)) }))
chk("Var*(IPW score) = Gamma/r (MC, relative)", max(abs(cov(sims) - G / 30)) / max(abs(G / 30)), 0.1)
## Poisson design: exact conditional variance and its HT estimator
pp <- poisson_inclusion(sqrt(rowSums(q^2)), st$ev, r = 150)$p
Vex <- crossprod(q * sqrt((1 - pp) / pp)) / st$n^2
sims <- t(replicate(20000, { R <- runif(st$n) < pp; colSums(q * (R / pp)) / st$n }))
chk("Poisson Var*(score) = m^-2 sum (1-p)/p q q^T", max(abs(cov(sims) - Vex)) / max(abs(Vex)), 0.1)
ht <- Reduce(`+`, replicate(4000, { R <- runif(st$n) < pp; crossprod(q * sqrt(R * (1 - pp) / pp^2)) / st$n^2 }, simplify = FALSE)) / 4000
chk("E*(HT variance estimator) = exact variance", max(abs(ht - Vex)) / max(abs(Vex)), 0.05)
chk("Poisson expected size = r", abs(sum(pp) - 150), 1e-8)
## case-cohort estimator: with all controls sampled it equals the full-cohort fit, and its
## estimate equals the explicit Kulich-Lin formula on a random case-cohort sample
fullfit <- ah_fit_full(st)
cc1 <- ah_fit_casecohort(st, rep(1, st$n), rep(0, st$n))
chk("case-cohort with p = 1 equals full fit", max(abs(cc1$beta - fullfit$beta)), 1e-10)
pc <- rep(0.3, st$n); Uc <- runif(st$n); rr <- ifelse(st$ev == 1, 1, (Uc < pc) / pc)
cc <- ah_fit_casecohort(st, pc, Uc)
tg <- st$ob$tg; Xs <- st$X; tt_ <- st$ob$t
Dm <- 0; nu <- 0
for (k in seq_along(tg)) {
  len <- min(tg[k], des$tau) - (if (k == 1) 0 else min(tg[k - 1], des$tau)); risk <- tt_ >= tg[k]
  if (sum(rr[risk]) <= 0) next
  zb <- colSums(Xs[risk, , drop = FALSE] * rr[risk]) / sum(rr[risk])
  if (len > 0) Dm <- Dm + len * crossprod(sweep(Xs[risk, , drop = FALSE], 2, zb) * sqrt(rr[risk]))
  evk <- which(tt_ == tg[k] & st$ev == 1); if (length(evk)) nu <- nu + colSums(sweep(Xs[evk, , drop = FALSE], 2, zb)) }
chk("case-cohort estimate equals explicit formula", max(abs(drop(solve(Dm, nu)) - cc$beta)), 1e-6)
## the GSELO screen does not depend on the unit of time; the split-full refit is equivariant
dsc <- calibrate_design(make_design(n = 40000, p = 20, cens = "uniform", cens_rate = 0.5), n_cal = 2e4)
dd <- gen_data(dsc, seed = 5)
g1 <- screen_gselo(dd$time, dd$status, dd$X, dsc$tau); g2 <- screen_gselo(dd$time * 30, dd$status, dd$X, dsc$tau * 30)
chk("GSELO selected set invariant to time unit", as.numeric(!identical(g1$A, g2$A)), 0.5)
g3 <- screen_gselo(dd$time, dd$status, dd$X, dsc$tau, ebic_n = "n")
g4 <- screen_gselo(dd$time * 30, dd$status, dd$X, dsc$tau * 30, ebic_n = "n")
chk("GSELO (log n criterion) invariant to time unit", as.numeric(!identical(g3$A, g4$A)), 0.5)
f1 <- ah_fit_full(ah_setup(dd$time, dd$status, dd$X[, dsc$S0], dsc$tau))
f2 <- ah_fit_full(ah_setup(dd$time * 30, dd$status, dd$X[, dsc$S0], dsc$tau * 30))
chk("split-full refit equivariant to time unit", max(abs(f1$beta - 30 * f2$beta) / f1$se), 1e-4)
## one split of the data-analysis protocol
sp <- primary_split(40000, 1)
a1 <- analyze_split(dd$time, dd$status, dd$X, dsc$tau, sp, "gselo", r_vec = c(500, 4000))
chk("analyze_split returns all methods", as.numeric(!all(c("SplitFull", "UIPW", "SIPW", "IIPW", "PUIPW", "PIIPW") %in% a1$table$method)), 0.5)
smoke_test()
if (!ok) stop("SOME CHECKS FAILED")
cat("ALL CHECKS PASSED\n")
