# How many subjects should BLAST keep (search_alignments)? When more subjects
# tie with the best one than BLAST returns (-max_target_seqs), the tied set is
# cut (`At best truncated`) and the LCA is computed on part of it.
# 30 seqs per marker from the benchmark pool, core_nt + the local 12S database,
# search_alignments 50 / 100 / 200, 30 cores; runs in the order
# 200 100 50 50 100 200 (warm cache; mean of 2 per setting).
# For each setting: wall time, seqs truncated, tied hits, and how the LCA
# compares with the one at 200. Writes only aggregated numbers.
# Run from repo root:
#   nohup nice -n 5 Rscript dev/benchmark/R/13_search_alignments.R > dev/benchmark/results/search_alignments.log 2>&1 &
suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(13)

db <- c(Sys.getenv("BLASTR_CORE_NT_DB"), Sys.getenv("BLASTR_LOCAL_DB"))
if (!all(nzchar(db))) stop("Set BLASTR_CORE_NT_DB and BLASTR_LOCAL_DB (see dev/benchmark/R/bench_lib.R).")
tax_db <- path.expand("~/.cache/R/BLASTr/taxonomy/taxonomy_2026-11-12.sqlite")
cores <- 30L
markers <- c("12S", "16S", "COI", "ITS2", "rbcL", "trnL")
n_per_marker <- 30L
out_dir <- "dev/benchmark/results/search_alignments"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pool <- read.delim("dev/benchmark/data/seqs_pool.tsv")
fa <- readLines("dev/benchmark/data/seqs_pool.fasta")
pool$sequence <- fa[match(paste0(">", pool$seq_id), fa) + 1L]
pool <- pool[pool$marker %in% markers & !is.na(pool$sequence), ]
pool <- do.call(rbind, lapply(split(pool, pool$marker), function(x) x[sample(nrow(x), min(nrow(x), n_per_marker)), ]))
pool <- pool[!duplicated(pool$sequence), ]
message(format(Sys.time()), " seqs: ", paste(names(table(pool$marker)), table(pool$marker), collapse = ", "))

runs <- c(200L, 100L, 50L, 50L, 100L, 200L)
res <- list()
times <- data.frame()
for (r in seq_along(runs)) {
  sa <- runs[r]
  t0 <- Sys.time()
  x <- parallel_blast(pool$sequence, db, total_cores = cores, search_alignments = sa, verbose = "silent")
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  times <- rbind(times, data.frame(run = r, search_alignments = sa, wall_s = wall))
  message(format(Sys.time()), " run ", r, " search_alignments ", sa, ": ", round(wall), " s")
  if (is.null(res[[as.character(sa)]])) res[[as.character(sa)]] <- x
}
saveRDS(list(res = res, times = times, marker = stats::setNames(pool$marker, pool$sequence)), file.path(out_dir, "runs.rds"))

per_seq <- do.call(rbind, lapply(names(res), function(sa) {
  x <- res[[sa]]
  s <- blast_summary(x)
  l <- blast_lca(x, taxonomy_db = tax_db)
  d <- dplyr::left_join(s[, c("Sequence", "N hits at best", "At best truncated")], l, by = "Sequence")
  d$search_alignments <- as.integer(sa)
  d
}))
per_seq$marker <- pool$marker[match(per_seq$Sequence, pool$sequence)]
ref <- per_seq[per_seq$search_alignments == 200L, c("Sequence", "LCA taxid", "LCA ranked rank", "LCA ranked name")]
names(ref)[-1] <- paste0(names(ref)[-1], " 200")
per_seq <- dplyr::left_join(per_seq, ref, by = "Sequence")

ranks <- c("species", "genus", "family", "order")
rank_level <- function(r) {
  i <- match(r, c("species", "subfamily", "genus", "family", "suborder", "order", "subclass", "class", "subphylum", "phylum", "kingdom", "domain"))
  i
}
agg <- function(d) data.frame(
  n = nrow(d),
  truncated = sum(d$`At best truncated`, na.rm = TRUE),
  median_at_best = stats::median(d$`N hits at best`, na.rm = TRUE),
  max_at_best = max(d$`N hits at best`, na.rm = TRUE),
  lca_species = sum(d$`LCA ranked rank` %in% "species"),
  lca_genus = sum(d$`LCA ranked rank` %in% "genus"),
  lca_family = sum(d$`LCA ranked rank` %in% c("subfamily", "family")),
  lca_above = sum(!is.na(d$`LCA ranked rank`) & !d$`LCA ranked rank` %in% c("species", "genus", "subfamily", "family")),
  lca_differs_from_200 = sum(!is.na(d$`LCA taxid`) & d$`LCA taxid` != d$`LCA taxid 200`, na.rm = TRUE),
  lca_too_deep_vs_200 = sum(rank_level(d$`LCA ranked rank`) < rank_level(d$`LCA ranked rank 200`), na.rm = TRUE)
)
by_m <- split(per_seq, list(per_seq$marker, per_seq$search_alignments), drop = TRUE)
sm <- do.call(rbind, lapply(by_m, function(d) cbind(marker = d$marker[1], search_alignments = d$search_alignments[1], agg(d))))
all <- do.call(rbind, lapply(split(per_seq, per_seq$search_alignments), function(d) cbind(marker = "all", search_alignments = d$search_alignments[1], agg(d))))
sm <- rbind(sm, all)
sm <- sm[order(sm$marker, sm$search_alignments), ]
tm <- stats::aggregate(wall_s ~ search_alignments, times, mean)
write.table(sm, file.path("dev/benchmark/results", "search_alignments.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(times, file.path("dev/benchmark/results", "search_alignments_times.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(times, row.names = FALSE)
print(tm, row.names = FALSE)
print(sm, row.names = FALSE)
message("Search alignments test done")
