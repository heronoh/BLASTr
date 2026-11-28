# Search for the calibration of the tie rules (15_tie_rules.R): the public
# references of 10_rank_calibration.R (GenBank sequences of known species,
# cut to the marker region) searched again in core_nt without their own
# accessions, now with 500 subjects (so the tied set is not cut, RE-02) and
# the fields needed to count differences in bases: nident, qstart, qend,
# qlen and btop (position by position, for the region two hits share), and
# stitle (voucher codes in the headers).
# Run from repo root, after 10_rank_calibration.R:
#   BLASTR_CORE_NT_DB=<core_nt prefix> nohup nice -n 5 Rscript dev/benchmark/R/14_tie_rules_search.R \
#     > dev/benchmark/results/tie_rules_search.log 2>&1 &
suppressMessages(devtools::load_all(".", quiet = TRUE))

db <- Sys.getenv("BLASTR_CORE_NT_DB")
if (!nzchar(db)) stop("Set BLASTR_CORE_NT_DB to the core_nt database (path and prefix).")
cores <- 30L
rc_dir <- "dev/benchmark/results/rank_calibration"
out_dir <- "dev/benchmark/results/tie_rules"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ref <- readRDS(file.path(rc_dir, "references.rds"))
neg <- file.path(rc_dir, "negative.bsl")
stopifnot(file.exists(neg))
hits_rds <- file.path(out_dir, "hits.rds")
if (!file.exists(hits_rds)) {
  t0 <- Sys.time()
  s <- parallel_blast_hsps(
    ref$query, db,
    total_cores = cores, max_target_seqs = 500L, max_hsps = 1L, perc_id = 70, perc_qcov_hsp = 50,
    fields = c(
      "saccver", "stitle", "pident", "length", "nident", "qstart", "qend", "qlen",
      "qcovhsp", "bitscore", "staxid", "btop"
    ),
    extra_args = c("-negative_seqidlist", neg),
    verbose = "silent"
  )
  saveRDS(s, hits_rds)
  message("search: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
}
s <- readRDS(hits_rds)
stopifnot(!any(s$saccver %in% ref$saccver))
message(format(Sys.time()), " hits: ", nrow(s), " for ", length(unique(s$Sequence)), " of ", nrow(ref), " queries")
