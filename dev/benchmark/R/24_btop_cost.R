# Cost of the tie rule by references in the long results (review of hub
# 2026-11-24-2335, item 5): rows, size in memory and on disk (rds) and time of
# blast_long() with and without btop, on the hits of the primer benchmark as
# parallel_blast() keeps them (identity >= 80, 100 subjects, 10 hits shown).
# Run from repo root: Rscript dev/benchmark/R/24_btop_cost.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(dplyr))
source("dev/benchmark/R/tie_rules_lib.R")

out_dir <- "dev/benchmark/results/primer_bench"
b <- bench_truth(readRDS(file.path(out_dir, "queries.rds")))
s <- readRDS(file.path(out_dir, "hits.rds"))
s <- s[s$Sequence %in% b$truth$seq, ]
num <- c("pident", "length", "nident", "qstart", "qend", "qcovhsp", "bitscore")
s[num] <- lapply(s[num], as.numeric)
s$qid <- match(s$Sequence, b$truth$seq)
s <- s[order(s$qid, s$hit), ]
s <- s[s$pident >= 80, ]
s <- s[stats::ave(seq_len(nrow(s)), s$qid, FUN = seq_along) <= 100L, ]
hits <- tibble::tibble(
  query = as.character(s$qid), subject = s$saccver, `subject header` = s$stitle, identity = s$pident,
  length = s$length, nident = s$nident, `query start` = s$qstart, `query end` = s$qend, qcovhsp = s$qcovhsp,
  bitscore = s$bitscore, staxid = trimws(sub(";.*", "", as.character(s$staxid))), btop = s$btop,
  blast_rank = as.integer(stats::ave(seq_len(nrow(s)), s$qid, FUN = seq_along)),
  `Similarity score` = similarity_score(s$pident, s$qcovhsp)
)
ids <- unique(hits$query)
seqs <- b$truth$seq[as.integer(ids)]
u <- uninformative_references()
run <- function(h) {
  t0 <- Sys.time()
  l <- blast_long(h, ids, seqs, num_alignments = 10L, uninformative = u)
  list(long = l, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
with_btop <- run(hits)
without <- run(dplyr::select(hits, -btop))
f1 <- tempfile(fileext = ".rds")
f2 <- tempfile(fileext = ".rds")
saveRDS(with_btop$long, f1)
saveRDS(dplyr::select(without$long, -at_best_refs, -btop), f2)
res <- data.frame(
  version = c("without btop (as before)", "with btop"),
  sequences = length(ids),
  hits_searched = nrow(hits),
  rows = c(nrow(without$long), nrow(with_btop$long)),
  memory_mb = round(c(utils::object.size(without$long), utils::object.size(with_btop$long)) / 2^20, 1),
  rds_mb = round(c(file.size(f2), file.size(f1)) / 2^20, 2),
  seconds = round(c(without$secs, with_btop$secs), 1)
)
print(res, row.names = FALSE)
write.table(res, "dev/benchmark/results/btop_cost.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
