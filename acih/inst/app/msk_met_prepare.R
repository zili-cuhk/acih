## =====================================================================
## app/msk_met_prepare.R -- MSK-MET (cBioPortal study msk_met_2021).
## Inputs in app/data/msk_met_2021/ (from github.com/cBioPortal/datahub,
## public/msk_met_2021; Git LFS files):
##   data_clinical_patient.txt, data_clinical_sample.txt,
##   data_mutations.txt (genomic analysis), data_gene_panel_matrix.txt,
##   and app/data/gene_panels/impact341.txt, impact410.txt, impact468.txt
##   (reference_data/gene_panels/data_gene_panel_impact*.txt).
## Outputs:
##   app/msk_met_clinical.rds  20 clinical candidates (baseline at sequencing)
##   app/msk_met_genomic.rds   the 20 clinical candidates plus nonsynonymous-
##                             mutation indicators for the genes assayed by
##                             all three IMPACT panels and mutated in >= 1%
##                             of patients
## Cohort: patients with overall-survival time > 0 (24,561 patients,
## 10,095 deaths, 10,025 deaths within 60 months).  Time in months.
## =====================================================================
read_cbio <- function(f) utils::read.delim(f, comment.char = "#", check.names = FALSE,
                                           stringsAsFactors = FALSE, na.strings = c("", "NA"))

## continuous candidates: optional log(1 + x), median imputation,
## standardisation, clipping at +-4 (the theory assumes bounded covariates)
zc <- function(x, log1p = FALSE) {
  x <- as.numeric(x); if (log1p) x <- log1p(pmax(x, 0))
  x[is.na(x)] <- stats::median(x, na.rm = TRUE)
  pmin(pmax((x - mean(x)) / stats::sd(x), -4), 4)
}

prepare_msk_met <- function(dir = "app/data/msk_met_2021", panel_dir = "app/data/gene_panels",
                            out_dir = "app", min_freq = 0.01) {
  pt <- read_cbio(file.path(dir, "data_clinical_patient.txt"))
  sa <- read_cbio(file.path(dir, "data_clinical_sample.txt"))
  sa <- sa[!duplicated(sa$PATIENT_ID), ]                       # one sample per patient
  df <- merge(pt, sa, by = "PATIENT_ID")
  df <- df[is.finite(df$OS_MONTHS) & df$OS_MONTHS > 0 & !is.na(df$OS_STATUS), ]
  ind <- function(x, lev) as.numeric(!is.na(x) & x %in% lev)
  cov <- data.frame(
    ## patient and tumour characteristics
    age_std            = zc(df$AGE_AT_SEQUENCING),
    male               = as.numeric(toupper(df$SEX) %in% "MALE"),
    race_white         = as.numeric(toupper(df$RACE) %in% "WHITE"),
    sample_metastasis  = ind(df$SAMPLE_TYPE, "Metastasis"),
    ## (metastasis and metastatic-site counts are not used: they summarize the
    ## clinical course and include metastases recorded after sequencing, the
    ## time origin of overall survival)
    ## genomic summaries
    fga_std            = zc(df$FGA),
    tmb_std            = zc(df$TMB_NONSYNONYMOUS, log1p = TRUE),
    msi_score_std      = zc(df$MSI_SCORE, log1p = TRUE),
    ## organ system (reference: Core GI and the remaining systems)
    organ_thoracic        = ind(df$ORGAN_SYSTEM, "Thoracic"),
    organ_genitourinary   = ind(df$ORGAN_SYSTEM, "Genitourinary"),
    organ_developmental_gi = ind(df$ORGAN_SYSTEM, "Developmental GI Tract"),
    organ_breast          = ind(df$ORGAN_SYSTEM, "Breast"),
    organ_gynecologic     = ind(df$ORGAN_SYSTEM, "Gynecologic"),
    organ_skin            = ind(df$ORGAN_SYSTEM, "Skin"),
    organ_head_neck       = ind(df$ORGAN_SYSTEM, "Head and Neck"),
    ## primary cancer type within organ system
    cancer_nsclc       = ind(df$CANCER_TYPE, "Non-Small Cell Lung Cancer"),
    cancer_colorectal  = ind(df$CANCER_TYPE, "Colorectal Cancer"),
    cancer_prostate    = ind(df$CANCER_TYPE, "Prostate Cancer"),
    cancer_bladder     = ind(df$CANCER_TYPE, "Bladder Cancer"),
    cancer_pancreatic  = ind(df$CANCER_TYPE, "Pancreatic Cancer"),
    cancer_endometrial = ind(df$CANCER_TYPE, "Endometrial Cancer"))
  status <- as.integer(substr(as.character(df$OS_STATUS), 1, 1))
  dat <- data.frame(id = df$PATIENT_ID, sample_id = df$SAMPLE_ID, time = as.numeric(df$OS_MONTHS),
                    status = status, cov, check.names = FALSE)
  stopifnot(qr(as.matrix(cov))$rank == ncol(cov))
  miss <- c(age = sum(is.na(df$AGE_AT_SEQUENCING)), fga = sum(is.na(df$FGA)),
            tmb = sum(is.na(df$TMB_NONSYNONYMOUS)), msi = sum(is.na(df$MSI_SCORE)))
  message("MSK-MET missing values (median-imputed): ", paste(names(miss), miss, sep = " = ", collapse = ", "))
  message(sprintf("MSK-MET: %d patients, %d deaths, %d deaths within 60 months, %d clinical candidates",
                  nrow(dat), sum(dat$status), sum(dat$status == 1 & dat$time <= 60), ncol(cov)))
  saveRDS(list(data = dat, covs = names(cov)), file.path(out_dir, "msk_met_clinical.rds"))

  mf <- file.path(dir, "data_mutations.txt")
  if (!file.exists(mf) || file.size(mf) < 1e6) {
    message("data_mutations.txt not found (or a Git LFS pointer): genomic analysis skipped")
    return(invisible(dat))
  }
  panel <- function(p) {
    l <- readLines(file.path(panel_dir, paste0(p, ".txt")), warn = FALSE)
    strsplit(sub("^gene_list:\\s*", "", grep("^gene_list:", l, value = TRUE)), "\t")[[1]]
  }
  common <- Reduce(intersect, lapply(c("impact341", "impact410", "impact468"), panel))
  mu <- utils::read.delim(mf, comment.char = "#", stringsAsFactors = FALSE, quote = "",
                          colClasses = "character")[, c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")]
  keep <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
            "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation")
  mu <- mu[mu$Variant_Classification %in% keep & mu$Tumor_Sample_Barcode %in% dat$sample_id &
           mu$Hugo_Symbol %in% common, ]
  mu$row <- match(mu$Tumor_Sample_Barcode, dat$sample_id)
  mu <- unique(mu[, c("row", "Hugo_Symbol")])
  frq <- table(mu$Hugo_Symbol) / nrow(dat)
  genes <- sort(names(frq)[frq >= min_freq])
  G <- matrix(0, nrow(dat), length(genes), dimnames = list(NULL, paste0("mut_", genes)))
  jj <- match(mu$Hugo_Symbol, genes); ok <- !is.na(jj)
  G[cbind(mu$row[ok], jj[ok])] <- 1
  datg <- cbind(dat, G)
  message(sprintf("MSK-MET genomic: %d genes on all IMPACT panels, %d mutated in >= %.0f%% of patients; p = %d",
                  length(common), length(genes), 100 * min_freq, ncol(cov) + length(genes)))
  saveRDS(list(data = datg, covs = c(names(cov), colnames(G)), mut_freq = frq[genes]),
          file.path(out_dir, "msk_met_genomic.rds"))
  invisible(datg)
}
if (sys.nframe() == 0L) prepare_msk_met()
