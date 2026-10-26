## =====================================================================
## File 12_figures.R -- publication figures (one visual system for the
## simulation study and the data analysis), base R graphics only.
##  * colour encodes the design family and lightness the allocation:
##      Poisson designs blue (P-OPT dark, P-UNIF light), with-replacement
##      informative designs orange (OPT dark, ZUO light), simple rules grey
##      (SNORM dark, UNIF light), case-cohort green (dashed);
##    every method also has its own plotting symbol and a direct label, so
##    no comparison relies on colour alone (palette checked for colour-
##    vision deficiency; within-family pairs differ in lightness);
##  * variance estimators: total sandwich (ink), inference-split only
##    (red), subsampling only (violet), with distinct symbols;
##  * white panels with a light grid and no axis lines, Monte Carlo bands
##    as shaded strips, bold panel letters, Helvetica throughout.
## =====================================================================

FIG <- list(
  ink = "#1a1a1a", ink2 = "#4d4d4b", muted = "#8f8e89", grid = "#e6e6e3", panel = "#ffffff",
  band = "#f6f6f4", mc = "#e2ebf6", ref = "#3d3d3b",
  meth = c(PIIPW = "#1f5ba8", PUIPW = "#5b93d3", IIPW = "#b03f15", ZIPW = "#d5801c",
           SIPW = "#4a4a4a", UIPW = "#8a8a8a", CaseCohort = "#17a07a", SplitFull = "#1a1a1a"),
  pch = c(PIIPW = 21, PUIPW = 22, IIPW = 24, ZIPW = 23, SIPW = 4, UIPW = 0, CaseCohort = 1, SplitFull = 3),
  lty = c(PIIPW = 1, PUIPW = 1, IIPW = 1, ZIPW = 1, SIPW = 1, UIPW = 1, CaseCohort = 2, SplitFull = 1),
  lwd = c(PIIPW = 2.3, PUIPW = 1.6, IIPW = 2.0, ZIPW = 1.6, SIPW = 1.5, UIPW = 1.5, CaseCohort = 1.6, SplitFull = 1.5),
  ## display labels (internal method identifiers are kept in the code and CSV files)
  lab = c(PIIPW = "P-OPT", IIPW = "OPT", PUIPW = "P-UNIF", ZIPW = "ZUO", CaseCohort = "KL-CC",
          SIPW = "SNORM", UIPW = "UNIF", SplitFull = "FULL", OracleFull = "ORACLE", Naive = "NAIVE", CIPW = "OPT-C"),
  var = c(total = "#1a1a1a", omega = "#d9423f", gamma = "#6250d6"),
  vpch = c(total = 21, omega = 1, gamma = 2),
  vlab = expression(total = "Total sandwich", omega = Omega * "-only", gamma = Gamma * "-only"),
  seq3 = c("#9cc3ec", "#3f7fcc", "#123d73"))

.fig_open <- function(file, width, height) {
  grDevices::pdf(file, width = width, height = height, family = "Helvetica", pointsize = 8.5, useDingbats = FALSE)
}
.fig_par <- function(mar = c(3.1, 3.4, 2.0, 0.8)) {
  graphics::par(mar = mar, mgp = c(1.9, 0.35, 0), tcl = 0, las = 1, cex.axis = 0.86, cex.lab = 0.92,
                col.axis = FIG$ink2, col.lab = FIG$ink, fg = FIG$muted, bty = "n", lend = 1, ljoin = 1)
}
## filled points: border in the series colour, fill white for open symbols
.pt <- function(x, y, m, cex = 0.95, col = FIG$meth[m], pch = FIG$pch[m], bg = NULL) {
  bg <- bg %||% (if (pch %in% 21:25) col else NA)
  graphics::points(x, y, pch = pch, col = if (pch %in% 21:25) "white" else col, bg = bg, cex = cex,
                   lwd = if (pch %in% 21:25) 0.6 else 1.2, xpd = NA)
}
## panel fill in user coordinates (log axes handled)
.bg <- function(col = FIG$panel) {
  u <- graphics::par("usr")
  x <- if (graphics::par("xlog")) 10^u[1:2] else u[1:2]; y <- if (graphics::par("ylog")) 10^u[3:4] else u[3:4]
  graphics::rect(x[1], y[1], x[2], y[2], col = col, border = NA)
}
.ux <- function() { u <- graphics::par("usr")[1:2]; if (graphics::par("xlog")) 10^u else u }
.uy <- function() { u <- graphics::par("usr")[3:4]; if (graphics::par("ylog")) 10^u else u }
.panel <- function(tag, title, at = NULL) {
  at <- at %||% .ux()[1]
  graphics::mtext(tag, side = 3, adj = 0, line = 0.55, font = 2, cex = 1.05, col = FIG$ink, at = at)
  graphics::mtext(title, side = 3, adj = 0, line = 0.55, cex = 0.9, col = FIG$ink,
                  at = at + graphics::strwidth(paste0(tag, "  "), cex = 1.05, font = 2))
}
.panel_log <- function(tag, title) {
  ## for log-scaled x axes compute the title offset in log units
  u <- graphics::par("usr")[1]
  graphics::mtext(tag, side = 3, adj = 0, line = 0.55, font = 2, cex = 1.05, col = FIG$ink, at = 10^u)
  w <- graphics::strwidth(paste0(tag, "  "), cex = 1.05, font = 2, units = "figure") /
    (graphics::par("plt")[2] - graphics::par("plt")[1]) * diff(graphics::par("usr")[1:2])
  graphics::mtext(title, side = 3, adj = 0, line = 0.55, cex = 0.9, col = FIG$ink, at = 10^(u + w))
}
.hgrid <- function(at) graphics::abline(h = at, col = FIG$grid, lwd = 0.6)
.vgrid <- function(at) graphics::abline(v = at, col = FIG$grid, lwd = 0.6)
.axis <- function(side, at, labels = at, ...) graphics::axis(side, at = at, labels = labels, lwd = 0, ...)
## spread direct labels vertically so that they do not overlap (in user units)
.spread <- function(y, gap, lo = -Inf, hi = Inf) {
  o <- order(y); z <- y[o]; n <- length(z)
  for (k in seq_len(n)[-1]) z[k] <- max(z[k], z[k - 1] + gap)
  shift <- (mean(y) - mean(z)); z <- z + min(0, shift)
  if (z[1] < lo) z <- z + (lo - z[1])
  if (z[n] > hi) { z[n] <- hi; for (k in rev(seq_len(n - 1))) z[k] <- min(z[k], z[k + 1] - gap) }
  z[order(o)]
}
## direct labels at line ends (log axes): leader line, coloured symbol, label in ink, in the right margin
.endlabels <- function(endx, endy, gap, cex = 0.78) {
  u <- graphics::par("usr"); xr <- 10^u[2]; lo <- u[3] * log(10) + 0.4 * gap; hi <- u[4] * log(10) - 0.4 * gap
  ly <- exp(.spread(log(endy), gap, lo, hi))
  m <- names(endy); xs <- xr * 10^(0.03 * diff(u[1:2])); xt <- xr * 10^(0.065 * diff(u[1:2]))
  graphics::segments(endx * 10^(0.012 * diff(u[1:2])), endy, xs / 10^(0.012 * diff(u[1:2])), ly,
                     col = grDevices::adjustcolor(FIG$meth[m], 0.8), lwd = 0.7, xpd = NA)
  for (k in seq_along(m)) .pt(xs, ly[k], m[k], cex = 0.95)
  graphics::text(xt, ly, FIG$lab[m], adj = c(0, 0.45), cex = cex, col = FIG$ink, xpd = NA)
}
.legend_methods <- function(methods, x = "bottom", inset = c(0, 0.005), cex = 0.82, lines = TRUE) {
  graphics::legend(x, legend = FIG$lab[methods], col = ifelse(FIG$pch[methods] %in% 21:25, "white", FIG$meth[methods]),
                   pt.bg = FIG$meth[methods], pch = FIG$pch[methods], pt.lwd = ifelse(FIG$pch[methods] %in% 21:25, 0.6, 1.2),
                   pt.cex = 1.15, horiz = TRUE, bty = "n", cex = cex, x.intersp = 0.7,
                   text.width = max(graphics::strwidth(FIG$lab[methods], cex = cex)) * 1.15,
                   text.col = FIG$ink, inset = inset, xpd = NA)
  if (lines) invisible(NULL)
}
.design_group <- function(d) {
  s <- sub("^S1_m\\d+k_", "", d)
  unname(c(U50 = "50% uniform censoring", CD50 = "50% covariate-dependent censoring", U90 = "8% events", RARE4 = "4% events",
           LIN = "Other settings, m = 25,000", S10 = "Other settings, m = 25,000", BIN = "Other settings, m = 25,000")[s])
}
.design_sub <- function(d) {
  s <- sub("^S1_m\\d+k_", "", d)
  other <- c(LIN = "Linear baseline hazard", S10 = "Ten active coefficients", BIN = "Binary covariates")
  ifelse(s %in% names(other), other[s], sprintf("m = %s,000", sub("^S1_m(\\d+)k_.*", "\\1", d)))
}

## ---------------------------------------------------------------------
## Calibration: grouped range-dot plot of ASE/ESE and coverage by setting.
## For each method, the bar spans the budgets and the point is their mean;
## the shaded strip is the 95% Monte Carlo band around the nominal value.
## ---------------------------------------------------------------------
fig_calibration2 <- function(avg, file, methods = c("UIPW", "ZIPW", "IIPW", "PIIPW", "CaseCohort"),
                             R = 1000, row_lab = NULL, width = 7, height = NULL, groups = NULL) {
  avg <- avg[avg$method %in% methods & !is.na(avg$r) & avg$n_cond >= 50, ]
  ds <- unique(avg$design); nd <- length(ds); nm <- length(methods)
  if (is.null(row_lab)) { row_lab <- .design_sub(ds); groups <- groups %||% .design_group(ds) }
  groups <- groups %||% rep("", nd)
  ## vertical positions: a header line for each group with more than one row
  rl <- rle(groups); pos <- numeric(nd); hdr <- list(); y <- 0; k <- 0
  for (g in seq_along(rl$values)) {
    multi <- nzchar(rl$values[g])
    if (multi) { y <- y + 0.85; hdr[[length(hdr) + 1]] <- list(y = y, lab = rl$values[g]) }
    for (j in seq_len(rl$lengths[g])) { k <- k + 1; y <- y + 1; pos[k] <- y
    }
    y <- y + 0.35
  }
  gstart <- cumsum(c(1, rl$lengths))[seq_along(rl$lengths)]
  ylim <- c(max(pos) + 0.75, 0.25)
  height <- height %||% (0.95 + 0.215 * max(pos))
  .fig_open(file, width, height); on.exit(grDevices::dev.off())
  graphics::par(oma = c(2.2, 0, 0, 0))
  graphics::layout(matrix(1:2, 1), widths = c(1.4, 1))
  off <- seq(-0.3, 0.3, length.out = nm); names(off) <- methods
  panel <- function(var, xlim, ref, band, tag, title, xlab, left) {
    .fig_par(mar = c(2.7, left, 2.0, 0.9))
    graphics::plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = "", yaxt = "n", xaxt = "n", xaxs = "i", yaxs = "i")
    .bg()
    ## alternate group strips
    for (g in seq_along(rl$values)) if (g %% 2 == 0) {
      i1 <- gstart[g]; i2 <- gstart[g] + rl$lengths[g] - 1
      top <- pos[i1] - 0.5 - (if (nzchar(rl$values[g])) 0.85 else 0); graphics::rect(xlim[1], top - 0.1, xlim[2], pos[i2] + 0.5 + 0.1,
                                                                                   col = FIG$band, border = NA)
    }
    graphics::rect(band[1], ylim[2], band[2], ylim[1], col = grDevices::adjustcolor(FIG$mc, 0.85), border = NA)
    at <- pretty(xlim, 5); at <- at[at > xlim[1] & at < xlim[2]]
    .vgrid(at); graphics::abline(v = ref, col = FIG$ref, lty = 2, lwd = 0.9)
    .axis(1, at = at, labels = formatC(at, format = "f", digits = 2))
    if (left > 2) {
      graphics::axis(2, at = pos, labels = row_lab, lwd = 0, cex.axis = 0.84, col.axis = FIG$ink, line = -0.2)
      for (h in hdr) graphics::mtext(h$lab, side = 2, at = h$y, line = left - 0.5, adj = 0, las = 1, cex = 0.72, font = 2,
                                    col = FIG$ink2)
    }
    for (m in methods) for (i in seq_len(nd)) {
      x <- avg[[var]][avg$design == ds[i] & avg$method == m]; if (!length(x)) next
      yy <- pos[i] + off[m]; xm <- mean(x)
      if (length(x) > 1) graphics::segments(max(min(x), xlim[1]), yy, min(max(x), xlim[2]), yy,
                                            col = grDevices::adjustcolor(FIG$meth[m], 0.75), lwd = 1.6)
      if (xm < xlim[1] || xm > xlim[2]) {
        xe <- if (xm < xlim[1]) xlim[1] else xlim[2]
        graphics::arrows(xe + (if (xm < xlim[1]) 1 else -1) * 0.012 * diff(xlim), yy, xe, yy, length = 0.04, col = FIG$meth[m], lwd = 1.2)
      } else .pt(xm, yy, m, cex = 1)
    }
    .panel(tag, title)
  }
  se <- 1.96 * sqrt(0.95 * 0.05 / R)
  panel("ratio", c(0.88, 1.12), 1, 1 + c(-1, 1) * 1.96 / sqrt(2 * R), "a", "ASE / ESE", "Ratio of estimated to empirical SE",
        if (any(nzchar(groups))) 8.2 else 10.5)
  panel("CP", c(0.915, 0.985), 0.95, 0.95 + c(-1, 1) * se, "b", "Coverage of nominal 95% intervals", "Coverage", 0.6)
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); graphics::plot.new()
  .legend_methods(methods)
}

## ---------------------------------------------------------------------
## Calibration summary (main text): one row per method and study; faint
## dots are the setting-by-budget cells, the bar spans them, the thick bar
## is their interquartile range and the large point their median; shaded
## strips are the 95% Monte Carlo bands of each study.
## studies: list of list(avg, methods, R, label)
## ---------------------------------------------------------------------
fig_calibration_summary <- function(studies, file, width = 7, height = NULL) {
  rows <- list(); y <- 0; hdr <- list()
  for (s in seq_along(studies)) {
    st <- studies[[s]]; y <- y + 0.9; hdr[[s]] <- list(y = y, lab = st$label, R = st$R, y0 = y + 0.5)
    for (m in st$methods) { y <- y + 1; rows[[length(rows) + 1]] <- list(s = s, m = m, y = y) }
    hdr[[s]]$y1 <- y + 0.5; y <- y + 0.4
  }
  ylim <- c(y + 0.3, 0.3); height <- height %||% (0.9 + 0.235 * y)
  .fig_open(file, width, height); on.exit(grDevices::dev.off())
  graphics::layout(matrix(1:2, 1), widths = c(1.32, 1))
  panel <- function(var, xlim, ref, bandfun, tag, title, xlab, left) {
    .fig_par(mar = c(2.8, left, 2.0, 0.9))
    graphics::plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = "", yaxt = "n", xaxt = "n", xaxs = "i", yaxs = "i")
    for (h in hdr) { b <- bandfun(h$R)
      graphics::rect(b[1], h$y0, b[2], h$y1, col = FIG$mc, border = NA) }
    at <- pretty(xlim, 5); at <- at[at > xlim[1] & at < xlim[2]]
    .vgrid(at); graphics::abline(v = ref, col = FIG$ref, lty = 2, lwd = 0.9)
    .axis(1, at = at, labels = formatC(at, format = "f", digits = 2))
    for (r in rows) {
      a <- studies[[r$s]]$avg; x <- a[[var]][a$method == r$m & !is.na(a$r) & a$n_cond >= 50]; x <- x[is.finite(x)]
      if (!length(x)) next
      xc <- pmin(pmax(x, xlim[1]), xlim[2]); q <- stats::quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
      jit <- (stats::runif(length(x)) - 0.5) * 0.42
      graphics::segments(max(min(x), xlim[1]), r$y, min(max(x), xlim[2]), r$y, col = grDevices::adjustcolor(FIG$meth[r$m], 0.55), lwd = 1)
      graphics::points(xc, r$y + jit, pch = 16, cex = 0.45, col = grDevices::adjustcolor(FIG$meth[r$m], 0.4))
      graphics::segments(q[1], r$y, q[3], r$y, col = FIG$meth[r$m], lwd = 4.2, lend = 1)
      .pt(q[2], r$y, r$m, cex = 1.45)
      if (any(x < xlim[1])) graphics::arrows(xlim[1] + 0.015 * diff(xlim), r$y, xlim[1], r$y, length = 0.04, col = FIG$meth[r$m], lwd = 1)
      if (any(x > xlim[2])) graphics::arrows(xlim[2] - 0.015 * diff(xlim), r$y, xlim[2], r$y, length = 0.04, col = FIG$meth[r$m], lwd = 1)
    }
    if (left > 2) {
      graphics::axis(2, at = vapply(rows, `[[`, 0, "y"), labels = FIG$lab[vapply(rows, `[[`, "", "m")], lwd = 0, cex.axis = 0.9,
                     col.axis = FIG$ink, line = -0.2)
      for (h in hdr) graphics::mtext(h$lab, side = 2, at = h$y, line = left - 0.5, adj = 0, las = 1, cex = 0.76, font = 2, col = FIG$ink2)
    }
    .panel(tag, title)
  }
  set.seed(20261026)
  panel("ratio", c(0.88, 1.12), 1, function(R) 1 + c(-1, 1) * 1.96 / sqrt(2 * R), "a", "ASE / ESE",
        "Ratio of estimated to empirical SE", 12.5)
  set.seed(20261026)
  panel("CP", c(0.915, 0.985), 0.95, function(R) 0.95 + c(-1, 1) * 1.96 * sqrt(0.95 * 0.05 / R), "b",
        "Coverage of nominal 95% intervals", "Coverage", 0.6)
}

## ---------------------------------------------------------------------
## Efficiency: interval length relative to the split-full refit against the
## sampling fraction r/m (log-log); faint points = settings, lines = medians.
## ---------------------------------------------------------------------
fig_efficiency2 <- function(avg, file, groups, titles, methods = c("UIPW", "ZIPW", "IIPW", "PUIPW", "PIIPW", "CaseCohort"),
                            width = 7, height = 5.6, ncol = 2) {
  .fig_open(file, width, height); on.exit(grDevices::dev.off())
  graphics::layout(matrix(seq_len(ncol * ceiling(length(groups) / ncol)), ncol = ncol, byrow = TRUE))
  for (g in seq_along(groups)) {
    x <- avg[avg$design %in% groups[[g]], ]
    sf <- stats::setNames(x$AIL[x$method == "SplitFull"], x$design[x$method == "SplitFull"])
    x <- x[x$method %in% methods & !is.na(x$r) & x$n_cond >= 50, ]
    x$m <- as.numeric(sub(".*_m(\\d+)k_.*", "\\1", x$design)) * 1000
    x$m[is.na(x$m)] <- as.numeric(sub("^S4_.*", "250000", x$design[is.na(x$m)]))
    x$f <- x$r / x$m; x$rel <- x$AIL / sf[x$design]
    xl <- range(x$f) * c(0.78, 1.18); yl <- c(0.86, max(x$rel) * 1.1)
    .fig_par(mar = c(3.1, 3.6, 2.0, 5.2))
    graphics::plot(NA, xlim = xl, ylim = yl, log = "xy", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i",
                   xlab = "Sampling fraction r / m", ylab = "Interval length relative to FULL")
    .bg()
    yat <- c(1, 1.25, 1.5, 2, 3, 4, 6, 8); yat <- yat[yat <= yl[2]]
    xat <- sort(unique(signif(x$f, 2)))
    .hgrid(yat); .vgrid(xat)
    .axis(1, at = xat, labels = paste0(100 * xat, "%")); .axis(2, at = yat, labels = as.character(yat))
    graphics::abline(h = 1, col = FIG$ref, lty = 2, lwd = 0.9)
    endy <- c(); endx <- c()
    for (m in rev(methods)) {
      y <- x[x$method == m, ]; if (!nrow(y)) next
      graphics::points(y$f, y$rel, pch = 16, col = grDevices::adjustcolor(FIG$meth[m], 0.28), cex = 0.55)
      md <- tapply(y$rel, signif(y$f, 2), stats::median); fx <- as.numeric(names(md))
      graphics::lines(fx, md, col = FIG$meth[m], lwd = FIG$lwd[m], lty = FIG$lty[m])
      .pt(fx, md, m, cex = 0.95)
      endx[m] <- max(fx); endy[m] <- md[length(md)]
    }
    .endlabels(endx, endy, 0.062 * diff(log(yl)))
    .panel_log(letters[g], titles[g])
  }
}

## ---------------------------------------------------------------------
## Failure mechanisms: (a) variance ablation, (b) screening challenges.
## ---------------------------------------------------------------------
fig_failure2 <- function(avg, file, method_B = "IIPW", extra_designs = c("S4_bb4", "S4_bb2"), width = 7, height = 3.3) {
  ab <- avg[avg$method %in% c("IIPW", "IIPW_OmegaOnly", "IIPW_GammaOnly") & !is.na(avg$r) & grepl("^S1_", avg$design), ]
  w <- merge(ab[ab$method == "IIPW", c("design", "r", "CP", "gamma_frac")],
             ab[ab$method == "IIPW_OmegaOnly", c("design", "r", "CP")], by = c("design", "r"), suffixes = c("", "_om"))
  gm <- ab[ab$method == "IIPW_GammaOnly", c("design", "r", "CP")]; names(gm)[3] <- "CP_gm"
  w <- merge(w, gm, by = c("design", "r"), all.x = TRUE)
  .fig_open(file, width, height); on.exit(grDevices::dev.off())
  graphics::par(oma = c(1.4, 0, 0, 0))
  graphics::layout(matrix(1:2, 1), widths = c(1, 1.08))
  .fig_par(mar = c(3.1, 3.6, 2.0, 4.8))
  graphics::plot(NA, xlim = c(0.18, 1.02), ylim = c(0.15, 1.03), xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i",
                 xlab = "Share of subsampling variance in the total", ylab = "Coverage of nominal 95% intervals")
  .bg(); .hgrid(seq(0.2, 1, 0.1)); .vgrid(seq(0.2, 1, 0.2))
  se <- 1.96 * sqrt(0.95 * 0.05 / 1000)
  graphics::rect(0.18, 0.95 - se, 1.02, 0.95 + se, col = grDevices::adjustcolor(FIG$mc, 0.9), border = NA)
  .axis(1, at = seq(0.2, 1, 0.2)); .axis(2, at = seq(0.2, 1, 0.2), labels = formatC(seq(0.2, 1, 0.2), format = "f", digits = 1))
  graphics::abline(h = 0.95, col = FIG$ref, lty = 2, lwd = 0.9)
  endy <- c()
  for (v in c("omega", "gamma", "total")) {
    y <- switch(v, total = w$CP, omega = w$CP_om, gamma = w$CP_gm); ok <- is.finite(y)
    if (sum(ok) > 4) { lf <- stats::lowess(w$gamma_frac[ok], y[ok], f = 0.7)
      graphics::lines(lf, col = grDevices::adjustcolor(FIG$var[v], 0.85), lwd = 1.8); endy[v] <- lf$y[length(lf$y)] }
    pc <- FIG$vpch[v]
    graphics::points(w$gamma_frac[ok], y[ok], pch = pc, col = if (pc %in% 21:25) "white" else FIG$var[v],
                     bg = FIG$var[v], cex = 0.8, lwd = if (pc %in% 21:25) 0.5 else 1)
  }
  ly <- .spread(endy, 0.07)
  vl <- list(total = "Total", omega = expression(Omega * "-only"), gamma = expression(Gamma * "-only"))
  for (v in names(endy)) graphics::text(1.045, ly[v], vl[[v]], adj = c(0, 0.5), col = FIG$var[v], cex = 0.84, font = 2, xpd = NA)
  graphics::text(0.21, 0.27, expression(Gamma * "-only: subsampling variance alone"), adj = c(0, 0), cex = 0.72, col = FIG$ink2)
  graphics::text(0.21, 0.205, expression(Omega * "-only: inference-split variance alone"), adj = c(0, 0), cex = 0.72, col = FIG$ink2)
  .panel("a", "Variance components (Study 1, OPT)")
  ## (b)
  s3 <- avg[(grepl("^S3_", avg$design) | avg$design %in% extra_designs) & avg$method == method_B & !is.na(avg$r) &
              !grepl("tv|frailty", avg$design), ]
  s3 <- s3[stats::ave(s3$r, s3$design, FUN = max) == s3$r, ]
  ord <- c("S3_block", "S3_t5", "S3_weak", "S3_bin", "S3_s10", "S3_U90", "S3_RARE4", extra_designs)
  s3 <- s3[match(intersect(ord, s3$design), s3$design), ]
  labs <- c(S3_block = "Block correlation", S3_t5 = "Bounded t5 margins", S3_weak = "Weak signals",
            S3_bin = "Binary covariates", S3_s10 = "Ten active coefficients", S3_U90 = "8% events",
            S3_RARE4 = "4% events", S4_bb4 = "4% events, n = 500,000", S4_bb2 = "2% events, n = 500,000")
  orc <- avg[avg$method == "OracleFull", ]; orc <- orc$CP[match(s3$design, orc$design)]
  ny <- nrow(s3); y <- seq_len(ny); y[y > 7] <- y[y > 7] + 0.4
  .fig_par(mar = c(3.1, 9.0, 2.0, 0.8))
  graphics::plot(NA, xlim = c(0, 1.03), ylim = c(max(y) + 0.6, 0.4), yaxt = "n", xaxt = "n", xaxs = "i", yaxs = "i",
                 xlab = "Proportion of replications", ylab = "")
  .bg(); .vgrid(seq(0, 1, 0.2))
  graphics::abline(h = 7.7, col = FIG$muted, lwd = 0.6, lty = 3)
  graphics::text(0.02, 7.95, "Study 4", adj = c(0, 0.5), cex = 0.68, col = FIG$ink2, font = 3)
  graphics::text(0.02, 0.62, "Study 3 (n = 100,000)", adj = c(0, 0.5), cex = 0.68, col = FIG$ink2, font = 3)
  graphics::abline(v = 0.95, col = FIG$ref, lty = 2, lwd = 0.9)
  .axis(1, at = seq(0, 1, 0.2), labels = formatC(seq(0, 1, 0.2), format = "f", digits = 1))
  graphics::axis(2, at = y, labels = labs[s3$design], lwd = 0, cex.axis = 0.84, col.axis = FIG$ink, line = -0.2)
  cpc <- ifelse(s3$n_cond >= 50, s3$CP, NA)
  graphics::rect(0, y - 0.22, s3$sure, y + 0.22, col = grDevices::adjustcolor(FIG$muted, 0.35), border = NA)
  graphics::segments(s3$CP_pipeline, y, ifelse(is.na(cpc), s3$CP_pipeline, cpc), y, col = FIG$muted, lwd = 1)
  graphics::points(orc, y, pch = 2, col = FIG$ink2, cex = 0.85, lwd = 1)
  graphics::points(s3$CP_pipeline, y, pch = 4, col = FIG$var["omega"], cex = 0.9, lwd = 1.5)
  graphics::points(cpc, y, pch = 21, col = "white", bg = FIG$meth["PIIPW"], cex = 1.05, lwd = 0.6)
  .panel("b", "Screening challenges (OPT, largest budget)")
  graphics::par(fig = c(0.47, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); graphics::plot.new()
  graphics::legend("bottom", legend = c("All active retained", "Conditional", "Pipeline", "ORACLE"),
                   pch = c(15, 21, 4, 2), col = c(grDevices::adjustcolor(FIG$muted, 0.5), "white", FIG$var["omega"], FIG$ink2),
                   pt.bg = c(NA, FIG$meth["PIIPW"], NA, NA), pt.cex = c(1.6, 1.1, 1, 0.9), ncol = 4,
                   bty = "n", cex = 0.74, text.col = FIG$ink, x.intersp = 0.5, inset = c(0, 0.005))
}

## ---------------------------------------------------------------------
## Data analysis: standard error relative to split-full against r.
## ---------------------------------------------------------------------
fig_app_efficiency2 <- function(eff_list, file, labels, methods = c("UIPW", "ZIPW", "SIPW", "IIPW", "PUIPW", "PIIPW", "CaseCohort"),
                                width = 7, height = 3.0) {
  .fig_open(file, width, height); on.exit(grDevices::dev.off())
  graphics::layout(matrix(seq_along(eff_list), 1))
  for (k in seq_along(eff_list)) {
    e <- eff_list[[k]]; e <- e[e$method %in% methods & !is.na(e$r), ]
    rr <- sort(unique(e$r)); yl <- c(0.86, max(e$se_rel) * 1.1)
    .fig_par(mar = c(3.1, if (k == 1) 3.6 else 2.4, 2.0, 4.6))
    graphics::plot(NA, xlim = range(rr) * c(0.8, 1.14), ylim = yl, log = "xy", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i",
                   xlab = "Subsample size r", ylab = if (k == 1) "Standard error relative to FULL" else "")
    .bg()
    yat <- c(1, 1.25, 1.5, 2, 3, 4, 6, 8); yat <- yat[yat <= yl[2]]
    .hgrid(yat); .vgrid(rr)
    .axis(1, at = rr, labels = paste0(rr / 1000, "k"), gap.axis = 0.1); .axis(2, at = yat, labels = as.character(yat))
    graphics::abline(h = 1, col = FIG$ref, lty = 2, lwd = 0.9)
    endy <- c(); endx <- c()
    for (m in rev(methods)) {
      y <- e[e$method == m, ]; if (!nrow(y)) next; y <- y[order(y$r), ]
      graphics::lines(y$r, y$se_rel, col = FIG$meth[m], lwd = FIG$lwd[m], lty = FIG$lty[m])
      .pt(y$r, y$se_rel, m, cex = 0.9)
      endx[m] <- max(y$r); endy[m] <- y$se_rel[nrow(y)]
    }
    .endlabels(endx, endy, 0.075 * diff(log(yl)), cex = 0.76)
    .panel_log(letters[k], labels[k])
  }
}

## ---------------------------------------------------------------------
## Plasmode: coverage of the true coefficients by variance estimator.
## ---------------------------------------------------------------------
fig_plasmode2 <- function(pav, file, labels, width = 7, height = 3.15, R = 500) {
  cohorts <- names(labels)
  .fig_open(file, width, height); on.exit(grDevices::dev.off())
  rows <- list(c("ZIPW", "ZIPW_GammaOnly", NA), c("IIPW", "IIPW_GammaOnly", "IIPW_OmegaOnly"), c("PIIPW", "PIIPW_GammaOnly", NA))
  rlab <- unname(FIG$lab[c("ZIPW", "IIPW", "PIIPW")])
  graphics::par(oma = c(1.6, 0, 0, 0))
  graphics::layout(matrix(seq_along(cohorts), 1), widths = c(1.42, rep(1, length(cohorts) - 1)))
  se <- 1.96 * sqrt(0.95 * 0.05 / R)
  for (k in seq_along(cohorts)) {
    x <- pav[pav$design == cohorts[k], ]; rr <- sort(unique(x$r[x$method == "PIIPW"]))
    nr <- length(rr); gap <- 0.7
    ypos <- lapply(seq_along(rows), function(g) (g - 1) * (nr + gap) + seq_len(nr))
    ymax <- max(unlist(ypos))
    .fig_par(mar = c(3.0, if (k == 1) 6.2 else 0.8, 2.0, 0.8))
    graphics::plot(NA, xlim = c(0, 1.02), ylim = c(ymax + 0.7, 0.3), yaxt = "n", xaxt = "n", xaxs = "i", yaxs = "i",
                   xlab = "Coverage of true coefficients", ylab = "")
    .bg()
    for (g in seq_along(rows)) if (g %% 2 == 0) graphics::rect(0, min(ypos[[g]]) - 0.5, 1.02, max(ypos[[g]]) + 0.5, col = FIG$band, border = NA)
    graphics::rect(0.95 - se, 0.3, 0.95 + se, ymax + 0.7, col = grDevices::adjustcolor(FIG$mc, 0.9), border = NA)
    .vgrid(seq(0, 1, 0.2)); graphics::abline(v = 0.95, col = FIG$ref, lty = 2, lwd = 0.9)
    .axis(1, at = seq(0, 1, 0.2), labels = formatC(seq(0, 1, 0.2), format = "f", digits = 1))
    for (g in seq_along(rows)) {
      for (j in seq_len(nr)) {
        i <- ypos[[g]][j]; r <- rr[j]
        for (v in c(3, 2, 1)) { m <- rows[[g]][v]; if (is.na(m)) next
          z <- x[x$method == m & x$r == r, ]; if (!nrow(z)) next
          vv <- c("total", "gamma", "omega")[v]
          graphics::segments(z$CP_min, i, z$CP_max, i, col = grDevices::adjustcolor(FIG$var[vv], 0.5), lwd = 1.5)
          pc <- FIG$vpch[vv]
          graphics::points(z$CP, i, pch = pc, col = if (pc %in% 21:25) "white" else FIG$var[vv], bg = FIG$var[vv],
                           cex = if (vv == "total") 1.4 else 0.95, lwd = if (pc %in% 21:25) 0.7 else 1.1) }
        if (k == 1) graphics::mtext(sprintf("r = %s", formatC(r, format = "d", big.mark = ",")), side = 2, at = i, line = 0.3,
                                    las = 1, adj = 1, cex = 0.68, col = FIG$ink2)
      }
      if (k == 1) graphics::mtext(rlab[g], side = 2, at = mean(ypos[[g]]), line = 4.4, las = 1, adj = 0.5, cex = 0.8,
                                  font = 2, col = FIG$ink)
    }
    .panel(letters[k], labels[k])
  }
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot.new()
  graphics::legend("bottom", legend = FIG$vlab, col = c("white", FIG$var[2:3]), pt.bg = FIG$var, pch = FIG$vpch,
                   pt.lwd = c(0.5, 1.1, 1.1), pt.cex = c(1.3, 1.1, 1.1), horiz = TRUE, bty = "n", cex = 0.9,
                   text.col = FIG$ink, x.intersp = 0.6, inset = c(0, 0.005))
}

## ---------------------------------------------------------------------
## Kaplan-Meier curves by score tertile, with numbers at risk.
## curves: data frame (cohort, group, time, surv); risk: (cohort, group, at, n)
## ---------------------------------------------------------------------
fig_km2 <- function(curves, risk, file, labels, units, width = 7, height = 3.5) {
  cohorts <- names(labels)
  .fig_open(file, width, height); on.exit(grDevices::dev.off())
  graphics::layout(matrix(seq_along(cohorts), 1), widths = c(1.32, rep(1, length(cohorts) - 1)))
  grp <- c("Low", "Middle", "High"); col <- stats::setNames(FIG$seq3, grp)
  for (k in seq_along(cohorts)) {
    cu <- curves[curves$cohort == cohorts[k], ]; ri <- risk[risk$cohort == cohorts[k], ]
    tau <- max(ri$at); ymin <- max(0, floor((min(cu$surv) - 0.02) * 20) / 20)
    .fig_par(mar = c(6.4, if (k == 1) 7.2 else 3.4, 2.0, 2.6))
    graphics::plot(NA, xlim = c(0, tau), ylim = c(ymin, 1.005), xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i",
                   xlab = sprintf("Time (%s)", units[k]), ylab = if (k == 1) "Survival probability" else "")
    .bg()
    xat <- pretty(c(0, tau), 4); xat <- xat[xat <= tau]; yat <- pretty(c(ymin, 1), 4); yat <- yat[yat >= ymin]
    .hgrid(yat); .vgrid(xat); .axis(1, at = xat); .axis(2, at = yat, labels = formatC(yat, format = "f", digits = 2))
    endy <- c()
    for (g in grp) { z <- cu[cu$group == g, ]
      graphics::lines(c(0, z$time), c(1, z$surv), type = "s", col = col[g], lwd = 2)
      endy[g] <- z$surv[nrow(z)] }
    ly <- .spread(endy, 0.06 * (1 - ymin))
    graphics::text(tau * 1.02, ly, grp, adj = c(0, 0.5), cex = 0.76, col = FIG$ink, xpd = NA)
    graphics::points(rep(tau * 1.005, 3), endy, pch = 21, col = "white", bg = col[grp], cex = 0.75, xpd = NA)
    for (j in seq_along(grp)) { g <- grp[j]; z <- ri[ri$group == g, ]
      graphics::mtext(formatC(z$n, format = "d", big.mark = ","), side = 1, line = 3.0 + 0.78 * j, at = z$at, cex = 0.6, col = FIG$ink2)
      if (k == 1) graphics::mtext(paste(g, "score"), side = 1, line = 3.0 + 0.78 * j, at = graphics::par("usr")[1] - 0.06 * tau,
                                  adj = 1, cex = 0.6, col = col[g], font = 2) }
    if (k == 1) graphics::mtext("Number at risk", side = 1, line = 3.0, at = graphics::par("usr")[1] - 0.06 * tau, adj = 1, cex = 0.6,
                                col = FIG$ink, font = 2)
    .panel(letters[k], labels[k])
  }
}

## ---------------------------------------------------------------------
## Redraw every figure from the summary CSV files (no refitting):
##   sim_dir: coefficient_averaged.csv; app_dir: <cohort>_efficiency.csv,
##   km_curves.csv, km_risk.csv, plasmode_averaged.csv.  The defaults follow
##   the layout of a run directory (results/ and app/results/).
## ---------------------------------------------------------------------
make_figures_from_csv <- function(dir = ".", sim_dir = file.path(dir, "results"), app_dir = file.path(dir, "app", "results"),
                                  out = file.path(dir, "results", "figures")) {
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(sim_dir, "coefficient_averaged.csv")
  if (file.exists(f)) make_sim_figures(utils::read.csv(f), out)
  lab <- c(msk_met = "MSK-MET clinical (p = 20)", msk_met_genomic = "MSK-MET genomic (p = 223)", nhis_lmf = "NHIS-LMF (3% mortality)")
  ef <- file.path(app_dir, paste0(names(lab), "_efficiency.csv")); ok <- file.exists(ef)
  if (any(ok)) {
    eff <- lapply(ef[ok], utils::read.csv); names(eff) <- names(lab)[ok]
    fig_app_efficiency2(eff, file.path(out, "fig_real_data_efficiency.pdf"), lab[ok])
  }
  kc <- file.path(app_dir, "km_curves.csv"); kr <- file.path(app_dir, "km_risk.csv")
  if (file.exists(kc) && file.exists(kr)) {
    cu <- utils::read.csv(kc); co <- intersect(names(lab), unique(cu$cohort))
    fig_km2(cu, utils::read.csv(kr), file.path(out, "fig_real_data_km.pdf"), lab[co],
            c(msk_met = "months", msk_met_genomic = "months", nhis_lmf = "years")[co])
  }
  pf <- file.path(app_dir, "plasmode_averaged.csv")
  if (file.exists(pf)) {
    pav <- utils::read.csv(pf); pl <- c(msk_met = "MSK-MET clinical", msk_met_genomic = "MSK-MET genomic", nhis_lmf = "NHIS-LMF")
    fig_plasmode2(pav, file.path(out, "fig_plasmode.pdf"), pl[intersect(names(pl), unique(pav$design))])
  }
  invisible(list.files(out, pattern = "^fig_.*\\.pdf$"))
}
