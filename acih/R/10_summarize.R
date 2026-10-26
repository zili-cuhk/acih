## 10_summarize.R -- tables (CSV + LaTeX) and figures (PDF) for Studies 1-3.

COLS <- c(SplitFull = "black", OracleFull = "grey50", Naive = "#a6761d", UIPW = "#1b9e77",
          SIPW = "#7570b3", ZIPW = "#66a61e", IIPW = "#d95f02", CIPW = "#e7298a",
          IIPW_OmegaOnly = "#e6ab02", RIPW = "#666666",
          PUIPW = "#1f78b4", PIIPW = "#b2182b", CaseCohort = "#6a3d9a")

read_design <- function(dir, study, id) {
  fs <- list.files(file.path(dir, "out", study, id), pattern = "^reps_.*rds$", full.names = TRUE)
  unlist(lapply(fs, function(f) readRDS(f)$reps), recursive = FALSE)
}

## coordinate-level summary: conditional on retention and finiteness, and
## pipeline coverage (miss or nonfinite = noncoverage)
## coordinate-level summary.  Conditional coverage is over replications that
## retained all active variables (sure screening) with finite intervals, as in
## the manuscript; pipeline coverage counts a screening miss or a nonfinite
## interval as noncoverage.  The case-cohort comparator is defined only when the
## event-complete design is feasible (number of events below r); replications
## where it is not are excluded and counted in `feasible`.
summ_coord <- function(reps, id) {
  ok <- Filter(function(r) is.null(r$error) && !isTRUE(r$gate_fail), reps)
  if (!length(ok)) return(NULL)
  D <- do.call(rbind, lapply(seq_along(ok), function(i) cbind(rep = i, ok[[i]]$res)))
  ## outputs of version 1.0.0 have no case-cohort rows when the selected set is
  ## empty: add them as screening misses (they are rows of PUIPW with no estimate)
  if (any(D$method == "CaseCohort")) {
    pu <- D[D$method == "PUIPW", ]; cc_key <- paste(D$rep, D$r)[D$method == "CaseCohort"]
    feas_r <- unique(D$r[D$method == "CaseCohort" & D$event_complete %in% TRUE])
    miss <- pu[!(paste(pu$rep, pu$r) %in% cc_key) & pu$r %in% feas_r, ]
    if (nrow(miss)) {
      miss$method <- "CaseCohort"; miss$event_complete <- NA
      miss$est <- NA_real_; miss$se <- NA_real_; miss$retained <- FALSE
      D <- rbind(D, miss)
    }
  }
  infeas <- D$method == "CaseCohort" & D$event_complete %in% FALSE
  D <- D[!infeas, ]
  if (!nrow(D)) return(NULL)
  grp <- paste(D$rep, D$method, D$r)
  D$allret <- stats::ave(D$retained, grp, FUN = all) == 1
  D$fin <- D$allret & is.finite(D$se) & is.finite(D$est)
  D$cov <- D$fin & abs(D$est - D$true) <= stats::qnorm(0.975) * D$se
  if (is.null(D$full)) D$full <- NA_real_
  D$cov_full <- D$fin & abs(D$est - D$full) <= stats::qnorm(0.975) * D$se   # coverage of the full-sample estimate
  key <- paste(D$method, D$r, D$coord)
  out <- do.call(rbind, lapply(split(D, key), function(x) {
    f <- x[x$fin, ]
    w <- wilson(sum(f$cov), max(nrow(f), 1))
    data.frame(design = id, method = x$method[1], r = x$r[1], coord = x$coord[1], true = x$true[1],
               R = nrow(x), feasible = nrow(x) / length(ok), sure = mean(x$allret),
               retained = mean(x$retained), n_cond = nrow(f), bias = mean(f$est - f$true),
               ESE = stats::sd(f$est), ASE = mean(f$se), ratio = mean(f$se) / stats::sd(f$est),
               CP = mean(f$cov), CP_lo = w[1], CP_hi = w[2], CP_pipeline = mean(x$cov),
               CP_full = if (all(is.na(f$full))) NA_real_ else mean(f$cov_full),
               AIL = mean(2 * stats::qnorm(0.975) * f$se), MSE = mean((f$est - f$true)^2),
               gamma_frac = mean(f$gamma_frac, na.rm = TRUE), ess = mean(f$ess, na.rm = TRUE),
               max_w = max(f$max_w, -Inf, na.rm = TRUE), time = mean(f$time, na.rm = TRUE))
  }))
  out[order(out$method, out$r, out$coord), ]
}
## coefficient-averaged summary (mean of coordinate-specific quantities)
summ_avg <- function(cs) {
  key <- paste(cs$design, cs$method, cs$r)
  do.call(rbind, lapply(split(cs, key), function(x)
    data.frame(design = x$design[1], method = x$method[1], r = x$r[1], R = x$R[1],
               feasible = x$feasible[1], sure = x$sure[1], n_cond = min(x$n_cond),
               retained = mean(x$retained), abs_bias = mean(abs(x$bias)), ESE = mean(x$ESE),
               ASE = mean(x$ASE), ratio = mean(x$ASE / x$ESE), CP = mean(x$CP),
               CP_min = min(x$CP), CP_max = max(x$CP), CP_pipeline = mean(x$CP_pipeline),
               CP_full = mean(x$CP_full),
               AIL = mean(x$AIL), MSE = mean(x$MSE), gamma_frac = mean(x$gamma_frac),
               time = mean(x$time))))
}
summ_screen <- function(reps, id) {
  ok <- Filter(function(r) is.null(r$error) && !isTRUE(r$gate_fail) && !is.null(r$screen), reps)
  if (!length(ok)) return(NULL)
  S <- do.call(rbind, lapply(ok, `[[`, "screen"))
  data.frame(design = id, R = nrow(S), size_mean = mean(S[, "size"]), size_min = min(S[, "size"]),
             size_max = max(S[, "size"]), sure = mean(S[, "sure"]), exact = mean(S[, "exact"]),
             boundary_low = mean(S[, "b_low"]), boundary_high = mean(S[, "b_high"]),
             converged = mean(S[, "conv"]), cap = mean(S[, "size"] >= 20), time = mean(S[, "t_screen"]),
             gate_fail = sum(sapply(reps, function(r) isTRUE(r$gate_fail))),
             errors = sum(sapply(reps, function(r) !is.null(r$error))))
}
summ_marginal <- function(reps, id) {
  ok <- Filter(function(r) is.null(r$error), reps)
  D <- do.call(rbind, lapply(ok, `[[`, "res"))
  do.call(rbind, lapply(split(D, D$K), function(x) {
    w <- wilson(sum(x$cand_sure), nrow(x)); m <- sum(x$cand_sure)
    data.frame(design = id, K = x$K[1], R = nrow(x), candidate = mean(x$cand_sure), cand_lo = w[1], cand_hi = w[2],
               gselo_cond = if (m) mean(x$gselo_sure[x$cand_sure], na.rm = TRUE) else NA, m = m,
               final = mean(x$cand_sure & x$gselo_sure %in% TRUE), gselo_fail = sum(x$gselo_fail),
               risk_gate = mean(x$risk_tau >= 0.05))
  }))
}

## ---- figures ----
## readable labels for design identifiers
design_label <- function(id) {
  lab <- c(U50 = "50% uniform cens.", CD50 = "50% cov.-dep. cens.", U90 = "90% cens. (8% events)",
           RARE4 = "no cens. (4% events)", LIN = "linear baseline", S10 = "10 active coef.",
           BIN = "binary covariates")
  parts <- strsplit(id, "_")
  vapply(parts, function(p) {
    if (p[1] == "S1") sprintf("m = %s, %s", sub("k", ",000", sub("^m", "", p[2])), lab[p[3]])
    else if (p[1] == "S2") sprintf("n = %s, p = %s%s", sub("k", ",000", sub("^n", "", p[2])), sub("^p", "", p[3]),
                                   if (p[4] == "CD50") ", cov.-dep." else "")
    else paste(p[-1], collapse = " ")
  }, "")
}
## budget level within a design: 1 = smallest r, ..., up to 3
budget_level <- function(avg) {
  lv <- stats::ave(avg$r, avg$design, FUN = function(r) match(r, sort(unique(r[!is.na(r)]))))
  lv[is.na(avg$r)] <- NA; lv
}
METH_COL <- c(UIPW = "#1b9e77", IIPW = "#d95f02", PUIPW = "#7570b3", PIIPW = "#e7298a",
              CaseCohort = "#666666", SplitFull = "black")

## Calibration (ASE/ESE and coverage) of selected methods by design and budget level
fig_calibration <- function(avg, file, methods = c("IIPW", "PIIPW"), R_ref = NULL) {
  a <- avg[avg$method %in% methods & !is.na(avg$r), ]
  a$lv <- budget_level(a)
  ids <- unique(a$design); y <- seq_along(ids)
  off <- stats::setNames(seq(-0.18, 0.18, length.out = length(methods)), methods)
  pch <- c(1, 2, 0)
  R_ref <- R_ref %||% max(a$R)
  grDevices::pdf(file, width = 8.5, height = 0.32 * length(ids) + 1.9)
  op <- graphics::par(mfrow = c(1, 2), mar = c(4, 12, 1.5, 0.8), oma = c(1.8, 0, 0, 0))
  for (v in c("ratio", "CP")) {
    xl <- if (v == "ratio") range(c(0.9, 1.1, a$ratio), na.rm = TRUE) else range(c(0.92, 0.98, a$CP), na.rm = TRUE)
    graphics::plot(NA, xlim = xl, ylim = rev(range(y)) + c(0.6, -0.6), yaxt = "n", ylab = "",
                   xlab = if (v == "ratio") "ASE / ESE" else "Coverage of 95% intervals")
    graphics::axis(2, at = y, labels = design_label(ids), las = 1, cex.axis = 0.72)
    graphics::abline(h = y, col = "grey92")
    if (v == "CP") { h <- 1.96 * sqrt(0.95 * 0.05 / R_ref)
      graphics::rect(0.95 - h, 0, 0.95 + h, length(ids) + 1, col = grDevices::adjustcolor("grey80", 0.5), border = NA) }
    graphics::abline(v = if (v == "ratio") 1 else 0.95, lty = 2)
    for (mth in methods) { b <- a[a$method == mth, ]
      graphics::points(b[[v]], match(b$design, ids) + off[mth], pch = pch[b$lv], col = METH_COL[mth], cex = 0.9) }
    graphics::mtext(if (v == "ratio") "(A)" else "(B)", side = 3, adj = 0, line = 0.3, font = 2)
  }
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  graphics::legend("bottom", c(methods, "smallest r", "middle r", "largest r"), xpd = TRUE, horiz = TRUE, bty = "n",
                   col = c(METH_COL[methods], rep("black", 3)), pch = c(rep(19, length(methods)), pch), cex = 0.75)
  graphics::par(op); grDevices::dev.off()
}

## Interval length relative to the split-full refit (cost of subsampling)
fig_efficiency <- function(avg, file, methods = c("UIPW", "IIPW", "PUIPW", "PIIPW", "CaseCohort")) {
  a <- avg[!is.na(avg$r) | avg$method == "SplitFull", ]
  sf <- a[a$method == "SplitFull", c("design", "AIL")]; names(sf)[2] <- "AIL_sf"
  a <- merge(a[a$method %in% methods, ], sf, by = "design")
  a$rel <- a$AIL / a$AIL_sf; a$lv <- budget_level(a)
  rare <- grepl("U90|RARE4", a$design)
  grDevices::pdf(file, width = 8.5, height = 3.9)
  op <- graphics::par(mfrow = c(1, 2), mar = c(4, 4.3, 1.8, 0.8), oma = c(1.8, 0, 0, 0))
  for (k in 1:2) {
    b <- a[if (k == 1) !rare else rare, ]
    yl <- range(c(1, b$rel), na.rm = TRUE)
    graphics::plot(NA, xlim = c(0.6, 3.4), ylim = yl, log = "y", xaxt = "n",
                   xlab = "Subsample size within a setting", ylab = "AIL / AIL(SplitFull)")
    graphics::axis(1, at = 1:3, labels = c("smallest r", "middle r", "largest r"))
    graphics::abline(h = 1, lty = 2)
    off <- stats::setNames(seq(-0.24, 0.24, length.out = length(methods)), methods)
    for (mth in methods) { cm <- b[b$method == mth & is.finite(b$rel), ]
      if (!nrow(cm)) next
      graphics::points(cm$lv + off[mth], cm$rel, col = grDevices::adjustcolor(METH_COL[mth], 0.45), pch = 19, cex = 0.6)
      md <- tapply(cm$rel, cm$lv, stats::median)
      graphics::lines(as.numeric(names(md)) + off[mth], md, col = METH_COL[mth], lwd = 2)
      graphics::points(as.numeric(names(md)) + off[mth], md, col = METH_COL[mth], pch = 19) }
    graphics::mtext(if (k == 1) "(A) 50% censoring settings" else "(B) Rare-event settings", side = 3, adj = 0, line = 0.3, font = 2)
  }
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  graphics::legend("bottom", methods, col = METH_COL[methods], lwd = 2, pch = 19, horiz = TRUE, bty = "n", cex = 0.8, xpd = TRUE)
  graphics::par(op); grDevices::dev.off()
}

## Failure mechanisms: (A) Omega-only versus total sandwich against the share of
## the subsampling component; (B) sure screening, conditional and pipeline coverage
## of IIPW (largest budget) in the screening challenges
fig_failure <- function(avg, file, method_B = "IIPW") {
  ab <- avg[avg$method %in% c("IIPW", "IIPW_OmegaOnly", "IIPW_GammaOnly") & !is.na(avg$r) & grepl("^S1_", avg$design), ]
  w <- merge(ab[ab$method == "IIPW", c("design", "r", "CP", "gamma_frac")],
             ab[ab$method == "IIPW_OmegaOnly", c("design", "r", "CP")], by = c("design", "r"), suffixes = c("", "_om"))
  gm <- ab[ab$method == "IIPW_GammaOnly", c("design", "r", "CP")]; names(gm)[3] <- "CP_gm"
  w <- merge(w, gm, by = c("design", "r"), all.x = TRUE)
  s3 <- avg[grepl("^S3_", avg$design) & avg$method == method_B & !is.na(avg$r) & !grepl("tv|frailty", avg$design), ]
  s3 <- s3[stats::ave(s3$r, s3$design, FUN = max) == s3$r, ]
  grDevices::pdf(file, width = 8.5, height = 3.9)
  op <- graphics::par(mfrow = c(1, 2), mar = c(4, 4.3, 1.8, 0.8))
  graphics::plot(w$gamma_frac, w$CP, ylim = c(0, 1), xlim = c(0, 1), pch = 19, col = METH_COL["IIPW"],
                 xlab = "Estimated share of the subsampling variance", ylab = "Coverage of 95% intervals")
  graphics::segments(w$gamma_frac, w$CP_om, w$gamma_frac, w$CP, col = "grey80")
  graphics::points(w$gamma_frac, w$CP_om, pch = 1)
  if (any(is.finite(w$CP_gm))) graphics::points(w$gamma_frac, w$CP_gm, pch = 2, col = "#1f78b4")
  graphics::abline(h = 0.95, lty = 2)
  graphics::legend("bottomleft", c("total sandwich", "inference-split variance only", "subsampling variance only"),
                   pch = c(19, 1, 2), col = c(METH_COL["IIPW"], "black", "#1f78b4"), bty = "n", cex = 0.75)
  graphics::mtext("(A) Variance ablation", side = 3, adj = 0, line = 0.3, font = 2)
  ids <- s3$design; y <- seq_along(ids)
  graphics::par(mar = c(4, 7.5, 1.8, 0.8))
  graphics::plot(NA, xlim = c(0, 1), ylim = rev(range(y)) + c(0.6, -0.6), yaxt = "n", ylab = "", xlab = "Proportion")
  graphics::axis(2, at = y, labels = c(block = "block corr.", t5 = "bounded t5", weak = "weak signals", s10 = "10 active",
                                        U90 = "8% events", RARE4 = "4% events", bin = "binary cov.")[sub("^S3_", "", ids)],
                 las = 1, cex.axis = 0.8)
  graphics::abline(h = y, col = "grey92"); graphics::abline(v = 0.95, lty = 2)
  orc <- avg[avg$method == "OracleFull", ]; orc <- orc$CP[match(ids, orc$design)]
  cpc <- ifelse(s3$n_cond >= 50, s3$CP, NA)          # conditional coverage needs >= 50 retained replications
  graphics::points(s3$sure, y, pch = 15, col = "grey40")
  graphics::points(cpc, y, pch = 19, col = METH_COL[method_B])
  graphics::points(s3$CP_pipeline, y, pch = 4, col = "black")
  graphics::points(orc, y, pch = 2, col = "grey30")
  graphics::legend("topleft", c("all active retained", "conditional coverage", "pipeline coverage", "OracleFull coverage"),
                   pch = c(15, 19, 4, 2), col = c("grey40", METH_COL[method_B], "black", "grey30"), bty = "n", cex = 0.72)
  graphics::mtext("(B) Screening challenges", side = 3, adj = 0, line = 0.3, font = 2)
  graphics::par(op); grDevices::dev.off()
}

latex_table <- function(df, cols, digits = 3, file) {
  f <- function(x) if (is.numeric(x)) ifelse(is.na(x), "--", formatC(x, digits = digits, format = "f")) else gsub("_", "\\\\_", x)
  lines <- c(sprintf("\\begin{tabular}{%s}", paste(rep("r", length(cols)), collapse = "")), "\\toprule",
             paste(gsub("_", "\\\\_", cols), collapse = " & "), "\\\\ \\midrule")
  lines[length(lines)] <- paste0(lines[3], " \\\\"); lines <- c(lines[1:2], lines[4], "\\midrule")
  for (i in seq_len(nrow(df))) lines <- c(lines, paste0(paste(sapply(cols, function(cn) f(df[[cn]][i])), collapse = " & "), " \\\\"))
  writeLines(c(lines, "\\bottomrule", "\\end{tabular}"), file)
}

summarize_results <- function(dir = ".") {
  RES <- file.path(dir, "results"); dir.create(RES, showWarnings = FALSE, recursive = TRUE)
  designs <- readRDS(file.path(dir, "out", "designs.rds"))
  coord <- list(); avg <- list(); scr <- list(); marg <- list()
  for (id in names(designs)) {
    des <- designs[[id]]; reps <- read_design(dir, des$study, id); if (!length(reps)) next
    if (des$run == "marginal") { marg[[id]] <- summ_marginal(reps, id); next }
    cs <- summ_coord(reps, id); if (is.null(cs)) { message("no usable replicates: ", id); next }
    coord[[id]] <- cs; avg[[id]] <- summ_avg(cs)
    if (des$run == "pipeline") scr[[id]] <- summ_screen(reps, id)
  }
  coord <- do.call(rbind, coord); avg <- do.call(rbind, avg)
  scr <- do.call(rbind, scr); marg <- do.call(rbind, marg)
  utils::write.csv(coord, file.path(RES, "coordinate_level.csv"), row.names = FALSE)
  utils::write.csv(avg, file.path(RES, "coefficient_averaged.csv"), row.names = FALSE)
  if (!is.null(scr)) utils::write.csv(scr, file.path(RES, "screening.csv"), row.names = FALSE)
  if (!is.null(marg)) utils::write.csv(marg, file.path(RES, "marginal_route.csv"), row.names = FALSE)
  s1 <- avg[grepl("^S1_", avg$design), ]; s2 <- avg[grepl("^S2_", avg$design), ]
  if (nrow(s1)) { make_sim_figures(avg, RES)
                  latex_table(s1[s1$method %in% c("SplitFull", "UIPW", "IIPW", "IIPW_OmegaOnly", "PUIPW", "PIIPW", "CaseCohort"), ],
                              c("design", "method", "r", "abs_bias", "ESE", "ASE", "ratio", "CP", "AIL", "gamma_frac"),
                              file = file.path(RES, "table_S1.tex")) }
  if (nrow(s2)) {
                  latex_table(s2, c("design", "method", "r", "ESE", "ASE", "ratio", "CP", "CP_pipeline", "AIL", "time"),
                              file = file.path(RES, "table_S2.tex")) }
  if (!is.null(scr)) latex_table(scr, c("design", "size_mean", "sure", "exact", "boundary_low", "boundary_high", "time"),
                                 file = file.path(RES, "table_screening.tex"))
  invisible(list(coord = coord, avg = avg, screening = scr, marginal = marg))
}

## Figures of Section 4 and Web Appendix B (12_figures.R)
make_sim_figures <- function(avg, RES) {
  ord1 <- c("S1_m10k_U50", "S1_m25k_U50", "S1_m50k_U50", "S1_m10k_CD50", "S1_m25k_CD50", "S1_m50k_CD50", "S1_m10k_U90",
            "S1_m25k_U90", "S1_m50k_U90", "S1_m25k_RARE4", "S1_m50k_RARE4", "S1_m25k_LIN", "S1_m25k_S10", "S1_m25k_BIN")
  s1 <- avg[avg$design %in% ord1, ]; s1 <- s1[order(match(s1$design, ord1)), ]
  if (nrow(s1)) {
    fig_calibration2(s1, file.path(RES, "fig_S1_calibration.pdf"))
    g <- list(ord1[!grepl("U90|RARE4", ord1)], ord1[grepl("U90", ord1)], ord1[grepl("RARE4", ord1)])
    ti <- c("50% censoring (Study 1)", "8% events (Study 1)", "4% events (Study 1)")
    if (any(avg$design == "S4_bb4")) { g <- c(g, list("S4_bb4")); ti <- c(ti, "4% events, n = 500,000 (Study 4)") }
    fig_efficiency2(avg, file.path(RES, "fig_S1_efficiency.pdf"), g, ti)
    fig_failure2(avg, file.path(RES, "fig_failure_mechanisms.pdf"))
  }
  o2 <- c("S2_n50k_p500_U50", "S2_n50k_p1000_U50", "S2_n100k_p500_U50", "S2_n100k_p1000_U50", "S2_n200k_p500_U50",
          "S2_n200k_p1000_U50", "S2_n100k_p500_CD50")
  s2 <- avg[avg$design %in% o2, ]; s2 <- s2[order(match(s2$design, o2)), ]
  if (nrow(s1) && nrow(s2))
    fig_calibration_summary(list(
      list(avg = s1[s1$method %in% c("UIPW", "ZIPW", "IIPW", "PIIPW", "CaseCohort"), ], methods = c("UIPW", "ZIPW", "IIPW", "PIIPW", "CaseCohort"),
           R = 1000, label = "Study 1: refit, 14 settings"),
      list(avg = s2, methods = c("UIPW", "ZIPW", "IIPW", "PUIPW", "PIIPW"), R = 500, label = "Study 2: complete pipeline, 7 designs")),
      file.path(RES, "fig_calibration_summary.pdf"))
  if (nrow(s2)) fig_calibration2(s2, file.path(RES, "fig_S2_calibration.pdf"), methods = c("UIPW", "ZIPW", "IIPW", "PUIPW", "PIIPW"),
    R = 500, row_lab = c("n = 50,000, p = 500", "n = 50,000, p = 1,000", "n = 100,000, p = 500", "n = 100,000, p = 1,000",
                         "n = 200,000, p = 500", "n = 200,000, p = 1,000", "n = 100,000, p = 500")[o2 %in% s2$design],
    groups = c(rep("50% uniform censoring", 6), "Covariate-dependent censoring")[o2 %in% s2$design])
  invisible(TRUE)
}
