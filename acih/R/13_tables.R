## =====================================================================
## File 13_tables.R -- the LaTeX tables of the manuscript and its
## supplement, computed from the summary CSV files (base R only).
## make_tables_from_csv() writes one file per table, <label>.tex, holding
## the tabular environment (header and body); captions are in the
## manuscript.  Tables printed over several pages in the supplement
## (Web Tables 7-8, 9-10 and 23-25) are written as one table here.
## =====================================================================

.f <- function(x, d) ifelse(is.na(x) | !is.finite(x), "--", formatC(x, format = "f", digits = d))
.f3 <- function(x) .f(x, 3); .f2 <- function(x) .f(x, 2)
.fint <- function(r) ifelse(is.na(r), "--", formatC(round(r), format = "d", big.mark = ","))
.row <- function(...) paste0(paste(..., sep = " & "), "\\\\")
.mname <- function(m) {
  base <- sub("_(GammaOnly|OmegaOnly)$", "", m); lab <- FIG$lab[base]; lab[is.na(lab)] <- base[is.na(lab)]
  suf <- ifelse(grepl("_GammaOnly$", m), ", \\(\\Gamma\\) only", ifelse(grepl("_OmegaOnly$", m), ", \\(\\Omega\\) only", ""))
  paste0(unname(lab), suf)
}
.tabular <- function(cols, head, body, file) {
  writeLines(c(sprintf("\\begin{tabular}{%s}", cols), "\\toprule", head, "\\midrule", body, "\\bottomrule", "\\end{tabular}"), file)
  invisible(body)
}
## drop a trailing \addlinespace
.trim <- function(rows) { while (length(rows) && rows[length(rows)] == "\\addlinespace") rows <- rows[-length(rows)]; rows }
.pick <- function(x, m, r) { y <- x[x$method == m & (if (is.na(r)) is.na(x$r) else (!is.na(x$r) & x$r == r)), , drop = FALSE]; if (nrow(y)) y[1, ] else NULL }

TAB_S1_ORDER <- c("S1_m10k_U50", "S1_m25k_U50", "S1_m50k_U50", "S1_m10k_CD50", "S1_m25k_CD50", "S1_m50k_CD50", "S1_m10k_U90",
                  "S1_m25k_U90", "S1_m50k_U90", "S1_m25k_RARE4", "S1_m50k_RARE4", "S1_m25k_LIN", "S1_m25k_S10", "S1_m25k_BIN")
TAB_S2_ORDER <- c("S2_n50k_p500_U50", "S2_n50k_p1000_U50", "S2_n100k_p500_U50", "S2_n100k_p1000_U50", "S2_n200k_p500_U50",
                  "S2_n200k_p1000_U50", "S2_n100k_p500_CD50")
.s1lab <- function(d) {
  p <- strsplit(d, "_")[[1]]
  c(m = sub("k", ",000", sub("m", "", p[2])),
    set = unname(c(U50 = "50\\% uniform", CD50 = "50\\% cov.-dep.", U90 = "90\\% uniform", RARE4 = "none (4\\% events)",
                   LIN = "linear baseline", S10 = "ten active", BIN = "binary covariates")[p[3]]))
}
.s2lab <- function(d) { p <- strsplit(d, "_")[[1]]; c(n = sub("k", ",000", sub("n", "", p[2])), p = sub("p", "", p[3])) }

## variable labels of the data analyses
TAB_VN <- c(age_std = "Age (SD)", male = "Male", race_white = "White race", sample_metastasis = "Metastatic sample",
  fga_std = "Fraction genome altered (SD)", tmb_std = "TMB (log, SD)", msi_score_std = "MSI score (log, SD)",
  organ_thoracic = "Thoracic organ system", organ_genitourinary = "Genitourinary organ system",
  organ_developmental_gi = "Developmental GI organ system", organ_breast = "Breast organ system",
  organ_gynecologic = "Gynecologic organ system", organ_skin = "Skin organ system", organ_head_neck = "Head and neck organ system",
  cancer_nsclc = "Non-small cell lung cancer", cancer_colorectal = "Colorectal cancer", cancer_prostate = "Prostate cancer",
  cancer_bladder = "Bladder cancer", cancer_pancreatic = "Pancreatic cancer", cancer_endometrial = "Endometrial cancer",
  bmi_std = "Body mass index (SD)", self_rated_health_std = "Self-rated health (SD, worse)", married = "Married",
  myocardial_infarction = "Myocardial infarction", stroke = "Stroke", copd = "COPD", cancer_history = "Cancer history",
  diabetes = "Diabetes", hypertension = "Hypertension", chd = "Coronary heart disease", angina = "Angina", emphysema = "Emphysema",
  asthma = "Asthma", arthritis = "Arthritis", kidney = "Weak or failing kidneys", liver = "Liver condition",
  current_smoker = "Current smoker", former_smoker = "Former smoker", hispanic = "Hispanic", black = "Black race",
  bachelor_degree = "Bachelor's degree or higher")
.vname <- function(v) {
  out <- unname(TAB_VN[v]); g <- is.na(out) & startsWith(v, "mut_")
  out[g] <- sprintf("\\textit{%s} mutation", substring(v[g], 5)); out[is.na(out)] <- gsub("_", "\\\\_", v[is.na(out)]); out
}
TAB_COHORT <- c(msk_met = "MSK-MET clinical", msk_met_genomic = "MSK-MET genomic", nhis_lmf = "NHIS-LMF")

## ---------------------------------------------------------------------
## Simulation tables
## ---------------------------------------------------------------------
tables_simulation <- function(sim_dir, out) {
  a <- utils::read.csv(file.path(sim_dir, "coefficient_averaged.csv"))
  cl <- utils::read.csv(file.path(sim_dir, "coordinate_level.csv"))
  scr <- utils::read.csv(file.path(sim_dir, "screening.csv"))
  mar <- utils::read.csv(file.path(sim_dir, "marginal_route.csv"))
  s1 <- a[startsWith(a$design, "S1_"), ]
  rs_of <- function(x) sort(unique(x$r[!is.na(x$r)]))
  s1_hd <- "\\(m\\) & Setting & \\(r\\) & FULL & UNIF & SNORM & ZUO & OPT & OPT-C & P-UNIF & P-OPT & KL-CC\\\\"
  meths <- c("UIPW", "SIPW", "ZIPW", "IIPW", "CIPW", "PUIPW", "PIIPW", "CaseCohort")
  s1table <- function(val, fmt, rel = FALSE) {
    rows <- character()
    for (d in TAB_S1_ORDER) {
      lb <- .s1lab(d); x <- s1[s1$design == d, ]; sf <- x[x$method == "SplitFull", ]
      for (k in seq_along(rs_of(x))) {
        r <- rs_of(x)[k]
        v <- vapply(meths, function(m) { y <- .pick(x, m, r); if (is.null(y)) NA_real_ else if (rel) y$AIL / sf$AIL else y[[val]] }, 0)
        sfv <- if (rel) sf$AIL else sf[[val]]
        rows <- c(rows, do.call(.row, as.list(c(if (k == 1) lb[["m"]] else "", if (k == 1) lb[["set"]] else "", .fint(r),
                                                if (k == 1) fmt(sfv) else "", fmt(v)))))
      }
      rows <- c(rows, "\\addlinespace")
    }
    .trim(rows)
  }
  .tabular("rlrrrrrrrrrr", s1_hd, s1table("CP", .f3), file.path(out, "tab_supp_s1_cp.tex"))
  .tabular("rlrrrrrrrrrr", s1_hd, s1table("ratio", .f3), file.path(out, "tab_supp_s1_ratio.tex"))
  .tabular("rlrrrrrrrrrr", s1_hd, s1table("AIL", .f2, rel = TRUE), file.path(out, "tab_supp_s1_ail.tex"))
  ## variance ablation, contrast targeting, Gamma-only intervals
  rows_ab <- rows_ct <- rows_gm <- character()
  for (d in TAB_S1_ORDER) {
    lb <- .s1lab(d); x <- s1[s1$design == d, ]; xc <- cl[cl$design == d, ]
    for (k in seq_along(rs_of(x))) {
      r <- rs_of(x)[k]; first <- c(if (k == 1) lb[["m"]] else "", if (k == 1) lb[["set"]] else "", .fint(r))
      ii <- .pick(x, "IIPW", r); oo <- .pick(x, "IIPW_OmegaOnly", r)
      rows_ab <- c(rows_ab, do.call(.row, as.list(c(first, .f3(c(ii$gamma_frac, ii$CP, oo$CP, oo$ratio))))))
      y <- xc[!is.na(xc$r) & xc$r == r, ]
      g <- function(m, sel, col) y[[col]][y$method == m & sel(y$coord)]
      i1 <- g("IIPW", function(j) j == 1, "AIL"); c1 <- g("CIPW", function(j) j == 1, "AIL")
      io <- mean(g("IIPW", function(j) j > 1, "AIL")); co <- mean(g("CIPW", function(j) j > 1, "AIL"))
      rows_ct <- c(rows_ct, do.call(.row, as.list(c(first, .f3(c(c1 / i1, co / io, g("CIPW", function(j) j == 1, "CP")))))))
      v <- function(m, col) { z <- .pick(x, m, r); if (is.null(z)) NA_real_ else z[[col]] }
      rows_gm <- c(rows_gm, do.call(.row, as.list(c(first, .f3(c(v("IIPW", "gamma_frac"), v("ZIPW_GammaOnly", "CP"), v("ZIPW_GammaOnly", "CP_full"),
        v("IIPW_GammaOnly", "CP"), v("PIIPW_GammaOnly", "CP"), v("ZIPW", "CP"), v("IIPW", "CP"), v("PIIPW", "CP")))))))
    }
    rows_ab <- c(rows_ab, "\\addlinespace"); rows_ct <- c(rows_ct, "\\addlinespace"); rows_gm <- c(rows_gm, "\\addlinespace")
  }
  .tabular("rlrrrrr", "\\(m\\) & Setting & \\(r\\) & \\(\\Gamma\\) share & CP total & CP \\(\\Omega\\)-only & Ratio \\(\\Omega\\)-only\\\\",
           .trim(rows_ab), file.path(out, "tab_supp_variance_ablation.tex"))
  .tabular("rlrrrr", "\\(m\\) & Setting & \\(r\\) & AIL ratio, coef. 1 & AIL ratio, coefs. 2--4 & CP, coef. 1\\\\",
           .trim(rows_ct), file.path(out, "tab_supp_contrast.tex"))
  .tabular("rlrrrrrrrrr", paste0("\\(m\\) & Setting & \\(r\\) & \\(\\Gamma\\) share & \\multicolumn{2}{c}{ZUO, \\(\\Gamma\\)-only} & OPT, \\(\\Gamma\\) & P-OPT, \\(\\Gamma\\) & ZUO & OPT & P-OPT\\\\\n",
                                   "\\cmidrule(lr){5-6}\n & & & & \\(\\beta_0\\) & full & only & only & total & total & total\\\\"),
           .trim(rows_gm), file.path(out, "tab_supp_gamma_only.tex"))
  ## Study 2: screening (main text) and complete pipeline
  rows <- "\\multicolumn{9}{l}{\\(\\rho=0.50\\); 50\\% uniform censoring}\\\\"
  for (d in TAB_S2_ORDER) {
    if (endsWith(d, "CD50")) rows <- c(rows, "\\addlinespace", "\\multicolumn{9}{l}{\\(\\rho=0.75\\); 50\\% covariate-dependent censoring}\\\\")
    x <- scr[scr$design == d, ][1, ]; lb <- .s2lab(d)
    rows <- c(rows, .row(lb[["n"]], lb[["p"]], sprintf("%s (%d--%d)", .f2(x$size_mean), as.integer(x$size_min), as.integer(x$size_max)),
                         .f3(x$sure), .f3(x$exact), .f3(x$boundary_low), .f3(x$boundary_high), .f3(x$cap), .f(x$time, 1)))
  }
  .tabular("rrrrrrrrr", "\\(n\\) & \\(p\\) & Mean \\(|\\widehat S|\\) (range) & Sure & Exact & Low & High & Cap & Time\\\\", rows,
           file.path(out, "tab_sim_design.tex"))
  rows <- character()
  for (d in TAB_S2_ORDER) {
    lb <- .s2lab(d); x <- a[a$design == d, ]; first <- TRUE
    lab <- paste0(lb[["n"]], ", ", lb[["p"]], if (endsWith(d, "CD50")) ", cov.-dep." else "")
    for (m in c("SplitFull", "OracleFull", "Naive", "UIPW", "SIPW", "ZIPW", "IIPW", "PUIPW", "PIIPW", "ZIPW_GammaOnly"))
      for (r in if (m %in% c("SplitFull", "OracleFull", "Naive")) NA else c(1000, 5000)) {
        y <- .pick(x, m, r); if (is.null(y)) next
        rows <- c(rows, .row(if (first) lab else "", .mname(m), .fint(r), .f(y$abs_bias, 4), .f(y$ESE, 4), .f(y$ASE, 4), .f3(y$ratio),
                             .f3(y$CP), .f3(y$CP_pipeline), .f(y$AIL, 4), .f(y$time, 1))); first <- FALSE
      }
    rows <- c(rows, "\\addlinespace")
  }
  .tabular("llrrrrrrrrr", "Design (\\(n\\), \\(p\\)) & Method & \\(r\\) & Bias & ESE & ASE & Ratio & CP & Pipeline & AIL & Time\\\\",
           .trim(rows), file.path(out, "tab_supp_s2.tex"))
  ## Study 3: screening challenges, misspecification, marginal pre-screen
  lab3 <- c(S3_block = "Block exchangeable, \\(\\rho=0.75\\)", S3_t5 = "Bounded \\(t_5\\), Weibull cens.", S3_weak = "Weak signals",
            S3_s10 = "Ten active coefficients", S3_U90 = "90\\% censoring (8\\% events)", S3_RARE4 = "No censoring (4\\% events)",
            S3_bin = "Binary covariates")
  rows <- character()
  for (d in names(lab3)) {
    s <- scr[scr$design == d, ][1, ]; x <- a[a$design == d, ]; first <- TRUE
    for (m in c("SplitFull", "OracleFull", "UIPW", "IIPW", "PIIPW", "CaseCohort"))
      for (r in if (m %in% c("SplitFull", "OracleFull")) NA else rs_of(x[x$method == m, ])) {
        y <- .pick(x, m, r)
        if (is.null(y) || is.na(y$CP_pipeline) || (m == "CaseCohort" && y$feasible == 0)) next
        ok <- y$n_cond >= 50
        rows <- c(rows, .row(if (first) lab3[[d]] else "", if (first) .f3(s$sure) else "", if (first) .f(s$size_mean, 1) else "",
                             .mname(m), .fint(r), as.character(as.integer(y$n_cond)), if (ok) .f3(y$ratio) else "--",
                             if (ok) .f3(y$CP) else "--", .f3(y$CP_pipeline), if (ok) .f3(y$AIL) else "--")); first <- FALSE
      }
    rows <- c(rows, "\\addlinespace")
  }
  .tabular("lrrlrrrrrr", "Design & Sure & Mean \\(|\\widehat S|\\) & Method & \\(r\\) & Ret. & Ratio & CP & Pipeline & AIL\\\\",
           .trim(rows), file.path(out, "tab_supp_screening_challenges.tex"))
  rows <- character()
  for (d in c(S3_tv = "Time-varying coefficient", S3_frailty = "Shared frailty (dependent censoring)")) {
    dd <- names(which(c(S3_tv = "Time-varying coefficient", S3_frailty = "Shared frailty (dependent censoring)") == d))
    x <- a[a$design == dd, ]; first <- TRUE
    for (m in c("SplitFull", "UIPW", "SIPW", "ZIPW", "IIPW", "CIPW", "IIPW_OmegaOnly", "PUIPW", "PIIPW"))
      for (r in if (m == "SplitFull") NA else c(1000, 5000)) {
        y <- .pick(x, m, r); if (is.null(y)) next
        nm <- if (m == "IIPW_OmegaOnly") "OPT, \\(\\Omega\\)-only" else .mname(m)
        rows <- c(rows, .row(if (first) d else "", nm, .fint(r), .f(y$abs_bias, 4), .f3(y$ratio), .f3(y$CP), .f3(y$CP_min),
                             .f3(y$CP_max), .f3(y$AIL))); first <- FALSE
      }
    rows <- c(rows, "\\addlinespace")
  }
  .tabular("llrrrrrrr", "Setting & Method & \\(r\\) & Bias & Ratio & CP & CP min & CP max & AIL\\\\", .trim(rows),
           file.path(out, "tab_supp_misspec.tex"))
  lab4 <- c(S3_pn_indep_alt = "Independent, alternating", S3_pn_ar50_reinf = "AR(1) \\(\\rho=0.50\\), reinforcing",
            S3_pn_ar85_alt = "AR(1) \\(\\rho=0.85\\), alternating", S3_pn_block_reinf = "Block \\(\\rho=0.75\\), reinforcing",
            S3_pn_t5_reinf = "Bounded \\(t_5\\) \\(\\rho=0.50\\), reinforcing")
  rows <- character()
  for (d in names(lab4)) for (i in which(mar$design == d)) {
    x <- mar[i, ]
    rows <- c(rows, .row(if (x$K == 25) lab4[[d]] else "", as.character(as.integer(x$K)),
                         sprintf("%s (%s, %s)", .f3(x$candidate), .f3(x$cand_lo), .f3(x$cand_hi)),
                         if (x$m == 0) "--" else sprintf("%s (%d)", .f3(x$gselo_cond), as.integer(x$m)), .f3(x$final), .f3(x$risk_gate)))
  }
  .tabular("lrrrrr", "Design & \\(K_n\\) & Candidate (95\\% CI) & GSELO (count) & Final & Risk\\\\", rows,
           file.path(out, "tab_supp_marginal_prescreen_pgeqn.tex"))
  ## Study 4
  rows <- character()
  for (d in c(S4_bb4 = "About 4\\% events", S4_bb2 = "About 2\\% events")) {
    dd <- if (grepl("4", d)) "S4_bb4" else "S4_bb2"; x <- a[a$design == dd, ]; sfa <- x$AIL[x$method == "SplitFull"][1]; first <- TRUE
    for (m in c("SplitFull", "OracleFull", "UIPW", "SIPW", "ZIPW", "IIPW", "PUIPW", "PIIPW", "CaseCohort", "ZIPW_GammaOnly", "PIIPW_GammaOnly"))
      for (r in if (m %in% c("SplitFull", "OracleFull")) NA else c(10000, 20000, 50000)) {
        y <- .pick(x, m, r); if (is.null(y)) next
        ok <- y$n_cond >= 50
        ail <- if (m == "OracleFull") "1.00" else if (ok) .f2(y$AIL / sfa) else "--"
        rows <- c(rows, .row(if (first) d else "", .mname(m), .fint(r), if (ok) .f3(y$ratio) else "--", if (ok) .f3(y$CP) else "--",
                             .f3(y$CP_pipeline), if (ok) .f3(y$CP_full) else "--", ail)); first <- FALSE
      }
    rows <- c(rows, "\\addlinespace")
  }
  .tabular("llrrrrrr", "Setting & Method & \\(r\\) & Ratio & CP & Pipeline & Full & Rel. AIL\\\\", .trim(rows), file.path(out, "tab_supp_s4.tex"))
  ## screening criterion comparison
  f <- file.path(sim_dir, "screencmp.csv")
  if (file.exists(f)) {
    sc <- utils::read.csv(f)
    labc <- c(S2_n50k_p500_U50 = "Study 2: \\(n=50{,}000\\), \\(p=500\\)", S2_n50k_p1000_U50 = "Study 2: \\(n=50{,}000\\), \\(p=1000\\)",
              S2_n100k_p500_U50 = "Study 2: \\(n=100{,}000\\), \\(p=500\\)", S2_n100k_p1000_U50 = "Study 2: \\(n=100{,}000\\), \\(p=1000\\)",
              S2_n200k_p500_U50 = "Study 2: \\(n=200{,}000\\), \\(p=500\\)", S2_n200k_p1000_U50 = "Study 2: \\(n=200{,}000\\), \\(p=1000\\)",
              S2_n100k_p500_CD50 = "Study 2: covariate-dependent censoring", S3_block = "Study 3: block exchangeable",
              S3_t5 = "Study 3: bounded \\(t_5\\)", S3_weak = "Study 3: weak signals", S3_bin = "Study 3: binary covariates",
              S3_s10 = "Study 3: ten active", S3_U90 = "Study 3: about 8\\% events", S3_RARE4 = "Study 3: about 4\\% events",
              S4_bb4 = "Study 4: about 4\\% events", S4_bb2 = "Study 4: about 2\\% events")
    rows <- character()
    for (d in names(labc)) {
      n0 <- sc[sc$design == d & sc$criterion == "log n", ]; e0 <- sc[sc$design == d & sc$criterion == "log events", ]
      if (!nrow(n0) || !nrow(e0)) next
      rows <- c(rows, .row(labc[[d]], formatC(e0$n_event, format = "f", digits = 0, big.mark = ","), .f3(n0$sure), .f3(e0$sure),
                           .f3(n0$exact), .f3(e0$exact), .f2(n0$size), .f2(e0$size)))
    }
    .tabular("lrrrrrrr", paste0("Design & Events & \\multicolumn{2}{c}{Sure} & \\multicolumn{2}{c}{Exact} & \\multicolumn{2}{c}{Size}\\\\\n",
                                "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\n",
                                " & & \\(\\log n\\) & \\(\\log N\\) & \\(\\log n\\) & \\(\\log N\\) & \\(\\log n\\) & \\(\\log N\\)\\\\"),
             rows, file.path(out, "tab_supp_ebic_compare.tex"))
  }
  invisible(TRUE)
}

## ---------------------------------------------------------------------
## Data-analysis and plasmode tables
## ---------------------------------------------------------------------
tables_applications <- function(app_dir, out) {
  rd <- function(f) utils::read.csv(file.path(app_dir, f))
  cohorts <- names(TAB_COHORT)[file.exists(file.path(app_dir, paste0(names(TAB_COHORT), "_primary.csv")))]
  for (co in cohorts) {
    p <- rd(paste0(co, "_primary.csv"))
    rows <- .row(.vname(p$var), .f2(p$frequency), .f2(1000 * p$est), .f2(1000 * p$se),
                 sprintf("(%s, %s)", .f2(1000 * p$lo), .f2(1000 * p$hi)))
    fn <- if (co == "msk_met_genomic") "tab_real_data_primary.tex" else sprintf("tab_supp_primary_%s.tex", co)
    .tabular("lrrrr", "Variable & Frequency & Estimate & SE & 95\\% CI\\\\", rows, file.path(out, fn))
  }
  ## efficiency
  M <- c("UIPW", "SIPW", "IIPW", "PUIPW", "PIIPW", "CaseCohort"); rows <- character()
  for (co in cohorts) {
    e <- rd(paste0(co, "_efficiency.csv"))
    for (k in seq_along(sort(unique(e$r[!is.na(e$r)])))) {
      r <- sort(unique(e$r[!is.na(e$r)]))[k]
      v <- vapply(M, function(m) { z <- e$se_rel[e$method == m & !is.na(e$r) & e$r == r]; if (length(z)) z[1] else NA_real_ }, 0)
      rows <- c(rows, do.call(.row, as.list(c(if (k == 1) TAB_COHORT[[co]] else "", .fint(r), .f2(v)))))
    }
    rows <- c(rows, "\\addlinespace")
  }
  .tabular("lrrrrrrr", "Cohort & \\(r\\) & UNIF & SNORM & OPT & P-UNIF & P-OPT & KL-CC\\\\", .trim(rows),
           file.path(out, "tab_real_data_efficiency.tex"))
  ## screening summary
  rows <- character()
  for (co in cohorts) {
    s <- rd(paste0(co, "_screening.csv"))[1, ]; fq <- rd(paste0(co, "_selection_frequency.csv"))
    rows <- c(rows, .row(TAB_COHORT[[co]], nrow(fq), sprintf("%s (%d--%d)", .f(s$size_mean, 1), as.integer(s$size_min), as.integer(s$size_max)),
                         .f2(s$cap), sum(fq$frequency >= 0.95), formatC(s$events_inf, format = "f", digits = 0, big.mark = ","),
                         .f2(s$t_screen), .f2(s$neg_hazard)))
  }
  .tabular("lrrrrrrr", "Cohort & \\(p\\) & Mean \\(|\\widehat S|\\) (range) & Cap & Stable & Events & Time & Neg.\\\\", rows,
           file.path(out, "tab_supp_real_screening.tex"))
  ## size-cap sensitivity (genomic analysis)
  f20 <- file.path(app_dir, "msk_met_genomic_cap20_primary.csv")
  if (file.exists(f20)) {
    p40 <- rd("msk_met_genomic_primary.csv"); f40 <- rd("msk_met_genomic_selection_frequency.csv")
    p20 <- utils::read.csv(f20); if ("method" %in% names(p20)) p20 <- p20[p20$method == "SplitFull", ]
    q20 <- rd("msk_met_genomic_cap20_frequency.csv"); fr20 <- stats::setNames(q20[[if ("frequency" %in% names(q20)) "frequency" else "freq"]], q20$var)
    fr40 <- stats::setNames(f40$frequency, f40$var)
    v <- unique(c(p40$var, p20$var)); cl <- v[!startsWith(v, "mut_")]; ge <- sort(v[startsWith(v, "mut_")], method = "radix")
    v <- c(cl[order(.vname(cl), method = "radix")], ge)
    est <- function(p, x, col) { z <- p[[col]][p$var == x]; if (length(z)) .f2(1000 * z[1]) else "--" }
    fq <- function(f, x) .f2(if (is.na(f[x])) 0 else f[[x]])
    rows <- vapply(v, function(x) .row(.vname(x), fq(fr40, x), est(p40, x, "est"), est(p40, x, "se"), fq(fr20, x),
                                       est(p20, x, "est"), est(p20, x, "se")), "")
    .tabular("lrrrrrr", "Variable & Freq. (40) & Estimate (40) & SE (40) & Freq. (20) & Estimate (20) & SE (20)\\\\", unname(rows),
             file.path(out, "tab_supp_cap_genomic.tex"))
  }
  ## fixed-budget sensitivity
  for (co in cohorts) {
    f <- file.path(app_dir, paste0(co, "_budget_sensitivity.csv")); if (!file.exists(f)) next
    b <- utils::read.csv(f); b <- b[order(is.na(b$GSELO), method = "radix"), ]
    rows <- .row(.vname(b$var), .f2(1000 * b$GSELO), .f2(1000 * b$K5), .f2(1000 * b$K10), .f2(1000 * b$K15))
    .tabular("lrrrr", "Variable & GSELO & \\(K=5\\) & \\(K=10\\) & \\(K=15\\)\\\\", rows, file.path(out, sprintf("tab_supp_budget_%s.tex", co)))
  }
  ## survey-weighted sensitivity
  f <- file.path(app_dir, "nhis_lmf_survey_sensitivity.csv")
  if (file.exists(f)) {
    s <- utils::read.csv(f)
    .tabular("lrrrr", "Variable & U estimate & U SE & W estimate & W SE\\\\",
             .row(.vname(s$var), .f2(1000 * s$U_est), .f2(1000 * s$U_se), .f2(1000 * s$W_est), .f2(1000 * s$W_se)),
             file.path(out, "tab_supp_nhis_survey_sensitivity.tex"))
  }
  ## plasmode
  f <- file.path(app_dir, "plasmode_averaged.csv")
  if (file.exists(f)) {
    p <- utils::read.csv(f); info <- rd("plasmode_designs.csv"); rows <- character()
    for (co in intersect(names(TAB_COHORT), unique(p$design))) {
      x <- p[p$design == co, ]; sfa <- x$AIL[x$method == "SplitFull"][1]; first <- TRUE; ii <- info[info$cohort == co, ]
      lab <- sprintf("%s (\\(|S_0|=%d\\), \\(|A|=%d\\))", TAB_COHORT[[co]], as.integer(ii$n_support), as.integer(ii$n_A))
      for (m in c("SplitFull", "UIPW", "SIPW", "ZIPW", "IIPW", "PUIPW", "PIIPW", "CaseCohort", "ZIPW_GammaOnly", "IIPW_GammaOnly",
                  "PIIPW_GammaOnly", "IIPW_OmegaOnly"))
        for (r in if (m == "SplitFull") NA else sort(unique(x$r[x$method == m & !is.na(x$r)]))) {
          y <- .pick(x, m, r); if (is.null(y)) next
          rows <- c(rows, .row(if (first) lab else "", .mname(m), .fint(r), .f3(y$ratio), .f3(y$CP),
                               sprintf("%s--%s", .f3(y$CP_min), .f3(y$CP_max)), .f3(y$CP_full), .f2(y$AIL / sfa))); first <- FALSE
        }
      rows <- c(rows, "\\addlinespace")
    }
    .tabular("llrrrrrr", "Cohort & Method & \\(r\\) & Ratio & CP & CP range & Full & Rel. AIL\\\\", .trim(rows),
             file.path(out, "tab_supp_plasmode.tex"))
  }
  invisible(TRUE)
}

## Write every table of the manuscript and supplement (file names follow the
## table labels, e.g. tab_supp_s1_cp.tex for Web Table 1).  The defaults
## follow the layout of a run directory (results/ and app/results/).
make_tables_from_csv <- function(dir = ".", sim_dir = file.path(dir, "results"), app_dir = file.path(dir, "app", "results"),
                                 out = file.path(dir, "results", "tables")) {
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  if (file.exists(file.path(sim_dir, "coefficient_averaged.csv"))) tables_simulation(sim_dir, out)
  if (dir.exists(app_dir)) tables_applications(app_dir, out)
  invisible(list.files(out, pattern = "^tab_.*\\.tex$"))
}
