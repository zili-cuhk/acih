## =====================================================================
## app/msk_met_fetch_mutations.R -- mutation calls of the cBioPortal study
## msk_met_2021 through the public cBioPortal web API, written in the
## three-column format read by msk_met_prepare.R:
##   Hugo_Symbol <TAB> Tumor_Sample_Barcode <TAB> Variant_Classification
## The DataHub file data_mutations.txt (Git LFS) contains the same calls
## and can be used instead.  Needs internet access, the curl command-line
## tool and the jsonlite package (used only here).
## Usage (from the run directory):
##   Rscript app/msk_met_fetch_mutations.R [app/data/msk_met_2021/data_mutations.txt]
## =====================================================================
.cbio_call <- function(url, body = NULL) {
  tf <- tempfile(fileext = ".json")
  args <- c("-sS", "-L", "--retry", "3", "-o", tf, "-w", "%{http_code}", "-H", shQuote("Accept: application/json"))
  if (!is.null(body)) args <- c(args, "-X", "POST", "-H", shQuote("Content-Type: application/json"), "--data", shQuote(body))
  code <- suppressWarnings(as.integer(system2("curl", c(args, shQuote(url)), stdout = TRUE)))
  txt <- if (file.exists(tf)) readChar(tf, file.size(tf), useBytes = TRUE) else ""
  list(code = code, txt = txt)
}

fetch_msk_met_mutations <- function(out = "app/data/msk_met_2021/data_mutations.txt",
                                    base = "https://www.cbioportal.org/api", page_size = 10000) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("install the jsonlite package first: install.packages(\"jsonlite\")")
  prof <- "msk_met_2021_mutations"; sl <- "msk_met_2021_all"
  ## two equivalent endpoints: GET with the sample list, or POST .../fetch
  get_page <- function(k, how) {
    q <- sprintf("projection=DETAILED&pageSize=%d&pageNumber=%d", page_size, k)
    if (how == "GET") .cbio_call(sprintf("%s/molecular-profiles/%s/mutations?sampleListId=%s&%s", base, prof, sl, q))
    else .cbio_call(sprintf("%s/molecular-profiles/%s/mutations/fetch?%s", base, prof, q), sprintf("{\"sampleListId\":\"%s\"}", sl))
  }
  for (how in c("GET", "POST")) {
    pages <- list(); k <- 0; ok <- TRUE
    repeat {
      r <- get_page(k, how)
      if (is.na(r$code) || r$code != 200) {
        message(sprintf("%s request failed (HTTP %s): %s", how, r$code, substr(r$txt, 1, 300))); ok <- FALSE; break
      }
      m <- jsonlite::fromJSON(r$txt, flatten = TRUE)
      if (!length(m) || !NROW(m)) break
      pages[[k + 1]] <- data.frame(Hugo_Symbol = m[["gene.hugoGeneSymbol"]], Tumor_Sample_Barcode = m[["sampleId"]],
                                   Variant_Classification = m[["mutationType"]], stringsAsFactors = FALSE)
      message(sprintf("page %d: %d calls", k + 1, NROW(m)))
      if (NROW(m) < page_size) break
      k <- k + 1
    }
    if (ok && length(pages)) break
  }
  if (!ok || !length(pages)) stop("could not download the mutation calls; see the messages above")
  out_df <- do.call(rbind, pages)
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(out_df, out, sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("%d mutation calls in %d samples written to %s", nrow(out_df),
                  length(unique(out_df$Tumor_Sample_Barcode)), out))
  invisible(out_df)
}
if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  fetch_msk_met_mutations(if (length(a)) a[1] else "app/data/msk_met_2021/data_mutations.txt")
}
