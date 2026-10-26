## =====================================================================
## app/nhis_lmf_prepare.R -- NHIS 2016-2018 sample adult files linked to the
## public-use NCHS Linked Mortality Files (follow-up through 2019).
## Inputs in app/data/nhis/: samadult_<y>.csv, person_<y>.csv and
## NHIS_<y>_MORT_2019_PUBLIC.dat for y = 2016, 2017, 2018.
## LMF fixed-width layout (2019 public-use NHIS files): PUBLICID 1-14,
## ELIGSTAT 15, MORTSTAT 16, UCOD_LEADING 17-19, DIABETES 20, HYPERTEN 21,
## DODQTR 22, DODYEAR 23-26, WGT_NEW 27-34, SA_WGT_NEW 35-42.
## The printed counts should be 83,676 adults, 2,739 deaths and 2,545 deaths
## within three years, with 156 strata and 1,726 PSUs.
## =====================================================================
read_lmf <- function(f) {
  x <- utils::read.fwf(f, widths = c(14, 1, 1, 3, 1, 1, 1, 4, 8, 8), colClasses = "character",
                       strip.white = TRUE)
  names(x) <- c("PUBLICID", "ELIGSTAT", "MORTSTAT", "UCOD_LEADING", "DIABETES", "HYPERTEN",
                "DODQTR", "DODYEAR", "WGT_NEW", "SA_WGT_NEW")
  x
}
num <- function(x) suppressWarnings(as.numeric(x))

## sample-adult variables used below; person-file variables added when absent
## design variables: PSTRAT and PPSU from 2016 on (STRAT_P and PSU_P before)
SA_VARS <- c("HHX", "FMX", "FPX", "INTV_QRT", "WTFA_SA", "PSTRAT", "PPSU", "AGE_P", "SEX", "BMI",
             "R_MARITL", "MIEV", "STREV", "COPDEV", "CANEV", "DIBEV1", "HYPEV", "CHDEV", "ANGEV",
             "EPHEV", "AASMEV", "ARTH1", "KIDWKYR", "LIVYR", "SMKSTAT2", "HISPAN_I", "MRACRPI2")
PE_VARS <- c("PHSTAT", "EDUC1", "R_MARITL", "AGE_P", "SEX")

prepare_nhis <- function(dir = "app/data/nhis", years = 2016:2018, out_dir = "app") {
  all <- list()
  for (y in years) {
    sa <- utils::read.csv(file.path(dir, sprintf("samadult_%d.csv", y)), stringsAsFactors = FALSE)
    names(sa) <- toupper(names(sa))
    miss <- setdiff(SA_VARS, names(sa))
    pf <- file.path(dir, sprintf("person_%d.csv", y))
    pe <- utils::read.csv(pf, stringsAsFactors = FALSE); names(pe) <- toupper(names(pe))
    add <- setdiff(intersect(PE_VARS, names(pe)), names(sa))
    if (length(add)) sa <- merge(sa, pe[, c("HHX", "FMX", "FPX", add)], by = c("HHX", "FMX", "FPX"), all.x = TRUE)
    still <- setdiff(c(SA_VARS, PE_VARS), names(sa))
    message(sprintf("%d: %d sample adults; from person file: %s; missing: %s", y, nrow(sa),
                    paste(add, collapse = " "), if (length(still)) paste(still, collapse = " ") else "none"))
    for (v in still) sa[[v]] <- NA
    sa$PUBLICID <- sprintf("%d%06d%02d%02d", y, as.integer(sa$HHX), as.integer(sa$FMX), as.integer(sa$FPX))
    lm <- read_lmf(file.path(dir, sprintf("NHIS_%d_MORT_2019_PUBLIC.dat", y)))
    d <- merge(sa[, c("PUBLICID", SA_VARS, setdiff(PE_VARS, SA_VARS))], lm, by = "PUBLICID")
    message(sprintf("%d: %d linked records, %d eligible", y, nrow(d), sum(d$ELIGSTAT == "1")))
    d <- d[d$ELIGSTAT == "1", ]
    ## follow-up in years from the interview-quarter midpoint to the
    ## death-quarter midpoint or 31 December 2019
    t_int <- y + (num(d$INTV_QRT) - 0.5) / 4
    dead <- d$MORTSTAT == "1"
    t_end <- rep(2020, nrow(d))
    t_end[dead] <- num(d$DODYEAR[dead]) + (num(d$DODQTR[dead]) - 0.5) / 4
    d$time <- pmax(t_end - t_int, 1 / 8); d$status <- as.integer(dead)
    d$year <- y
    all[[as.character(y)]] <- d
  }
  common <- Reduce(intersect, lapply(all, names))
  d <- do.call(rbind, lapply(all, function(x) x[, common]))
  yn <- function(v) { v <- num(v); ifelse(v == 1, 1, ifelse(v == 2, 0, NA)) }
  z <- function(x) { x <- num(x); x[is.na(x)] <- stats::median(x, na.rm = TRUE); (x - mean(x)) / stats::sd(x) }
  mode_imp <- function(x) { if (all(is.na(x))) stop("variable entirely missing"); x[is.na(x)] <- as.numeric(names(which.max(table(x)))); x }
  bmi <- num(d$BMI) / 100; bmi[bmi > 99] <- NA
  phs <- num(d$PHSTAT); phs[phs > 5] <- NA
  educ <- num(d$EDUC1); educ[educ > 21] <- NA
  smk <- num(d$SMKSTAT2); smk[smk > 4] <- NA
  cov <- data.frame(
    age_std = z(d$AGE_P), male = as.numeric(num(d$SEX) == 1), bmi_std = z(bmi), self_rated_health_std = z(phs),
    married = mode_imp(ifelse(num(d$R_MARITL) %in% c(9), NA, as.numeric(num(d$R_MARITL) %in% c(1, 2, 3)))),
    myocardial_infarction = mode_imp(yn(d$MIEV)), stroke = mode_imp(yn(d$STREV)),
    copd = mode_imp(yn(d$COPDEV)), cancer_history = mode_imp(yn(d$CANEV)),
    diabetes = mode_imp(yn(d$DIBEV1)), hypertension = mode_imp(yn(d$HYPEV)),
    chd = mode_imp(yn(d$CHDEV)), angina = mode_imp(yn(d$ANGEV)), emphysema = mode_imp(yn(d$EPHEV)),
    asthma = mode_imp(yn(d$AASMEV)), arthritis = mode_imp(yn(d$ARTH1)),
    kidney = mode_imp(yn(d$KIDWKYR)), liver = mode_imp(yn(d$LIVYR)),
    current_smoker = mode_imp(ifelse(is.na(smk), NA, as.numeric(smk %in% c(1, 2)))),
    former_smoker = mode_imp(ifelse(is.na(smk), NA, as.numeric(smk == 3))),
    hispanic = mode_imp(as.numeric(num(d$HISPAN_I) != 12)),
    black = mode_imp(as.numeric(num(d$MRACRPI2) == 2)),
    bachelor_degree = mode_imp(ifelse(is.na(educ), NA, as.numeric(educ >= 18))))
  dat <- data.frame(time = d$time, status = d$status, cov)
  survey <- list(weight = num(d$WTFA_SA) / length(years), strata = paste(d$year, d$PSTRAT),
                 psu = paste(d$year, d$PSTRAT, d$PPSU))
  if (anyNA(d$PSTRAT) || anyNA(d$PPSU)) stop("design variables PSTRAT / PPSU missing")
  message(sprintf("NHIS-LMF: %d adults, %d deaths, %d deaths within 3 years, %d candidates; %d strata, %d PSUs",
                  nrow(dat), sum(dat$status), sum(dat$status == 1 & dat$time <= 3), ncol(cov),
                  length(unique(survey$strata)), length(unique(survey$psu))))
  message("candidate means: ", paste(sprintf("%s %.3f", names(cov), colMeans(cov)), collapse = "; "))
  saveRDS(list(data = dat, covs = names(cov), survey = survey), file.path(out_dir, "nhis_lmf.rds"))
  invisible(dat)
}
if (sys.nframe() == 0L) prepare_nhis()
