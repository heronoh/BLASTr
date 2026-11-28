# Multi-marker benchmark, step 2: search the queries of 17_primer_bench_queries.R
# in core_nt without the template entries (negative.bsl), as 14_tie_rules_search.R
# (500 subjects, btop). Distinct sequences are searched once.
# Run from repo root:
#   BLASTR_CORE_NT_DB=<core_nt prefix> nohup nice -n 5 Rscript dev/benchmark/R/18_primer_bench_search.R \
#     > dev/benchmark/results/primer_bench_search.log 2>&1 &
suppressMessages(devtools::load_all(".", quiet = TRUE))

db <- Sys.getenv("BLASTR_CORE_NT_DB")
if (!nzchar(db)) stop("Set BLASTR_CORE_NT_DB to the core_nt database (path and prefix).")
cores <- 30L
out_dir <- "dev/benchmark/results/primer_bench"
q <- readRDS(file.path(out_dir, "queries.rds"))
neg <- file.path(out_dir, "negative.bsl")
stopifnot(file.exists(neg))
hits_rds <- file.path(out_dir, "hits.rds")
if (!file.exists(hits_rds)) {
  t0 <- Sys.time()
  s <- parallel_blast_hsps(
    unique(q$seq), db,
    total_cores = cores, max_target_seqs = 500L, max_hsps = 1L, perc_id = 70, perc_qcov_hsp = 50,
    fields = c(
      "saccver", "stitle", "pident", "length", "nident", "qstart", "qend", "qlen",
      "qcovhsp", "bitscore", "staxid", "btop"
    ),
    extra_args = c("-negative_seqidlist", neg),
    verbose = "silent"
  )
  saveRDS(s, hits_rds)
  message("search: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}
s <- readRDS(hits_rds)
message(format(Sys.time()), " hits: ", nrow(s), " for ", length(unique(s$Sequence)), " of ", length(unique(q$seq)), " sequences")
