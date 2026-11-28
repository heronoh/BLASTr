# Multi-marker benchmark: checks that the package's tie rule by references
# (at_best_refs of parallel_blast(), blast_lca(ties = "at_best_refs") and its
# species flag) gives the results of the benchmark's obitag_floor_ident and
# obitag_floor_ident+gap (tie_rules_lib.R, 21_primer_bench_extra_rules.R) on
# the same hits, with the species in the database ("in") and out of it.
# Run from repo root: Rscript dev/benchmark/R/23_package_rule_check.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(dplyr))
source("dev/benchmark/R/tie_rules_lib.R")

out_dir <- "dev/benchmark/results/primer_bench"
tax_db <- path.expand("~/.cache/R/BLASTr/taxonomy/taxonomy_2026-11-12.sqlite") # taxonomy of core_nt_nov26
b <- bench_truth(readRDS(file.path(out_dir, "queries.rds")))
extra <- readRDS(file.path(out_dir, "rules_res_extra.rds"))

s <- readRDS(file.path(out_dir, "hits.rds"))
s <- s[s$Sequence %in% b$truth$seq, ]
num <- c("pident", "length", "nident", "qstart", "qend", "qlen", "qcovhsp", "bitscore")
s[num] <- lapply(s[num], as.numeric)
s$qid <- match(s$Sequence, b$truth$seq)
s$staxid <- trimws(sub(";.*", "", as.character(s$staxid)))
tax <- suppressWarnings(get_tax(unique(s$staxid), taxonomy_db = tax_db))
s$species <- tax$species[match(s$staxid, tax$taxid)]
s$true_species <- b$truth$species[s$qid]
# the hits in BLAST order, as default_setting(): identity >= 80, 100 subjects
s <- s[order(s$qid, s$hit), ]
s <- s[s$pident >= 80, ]
s$blast_rank <- as.integer(stats::ave(seq_len(nrow(s)), s$qid, FUN = seq_along))

check <- function(h, sc) {
  h <- h[stats::ave(seq_len(nrow(h)), h$qid, FUN = seq_along) <= 100L, ]
  h$blast_rank <- as.integer(stats::ave(seq_len(nrow(h)), h$qid, FUN = seq_along))
  hits <- tibble::tibble(
    query = as.character(h$qid), subject = h$saccver, `subject header` = h$stitle, identity = h$pident,
    length = h$length, nident = h$nident, `query start` = h$qstart, `query end` = h$qend, qcovhsp = h$qcovhsp,
    bitscore = h$bitscore, staxid = h$staxid, btop = h$btop, blast_rank = h$blast_rank,
    `Similarity score` = similarity_score(h$pident, h$qcovhsp)
  )
  ids <- unique(hits$query)
  long <- blast_long(hits, ids, b$truth$seq[as.integer(ids)], num_alignments = 10L, uninformative = uninformative_references())
  t0 <- Sys.time()
  lca <- blast_lca(long, ties = "at_best_refs", taxonomy = tax)
  message(sc, ": blast_lca ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s for ", length(ids), " queries")
  lca$qid <- match(lca$Sequence, b$truth$seq)
  r <- extra$runs[[sc]]
  ref <- r[r$rule == "obitag_floor_ident", c("qid", "rank", "species", "genus", "n_tied")]
  m <- merge(ref, lca[, c("qid", "LCA ranked rank", "LCA ranked name", "N tied", "Species flag")], by = "qid")
  same_call <- (m$rank %in% "species" & m$`LCA ranked rank` %in% "species" & m$species == m$`LCA ranked name`) |
    (m$rank %in% "genus" & m$`LCA ranked rank` %in% "genus" & m$genus == m$`LCA ranked name`) |
    (!m$rank %in% c("species", "genus") & !m$`LCA ranked rank` %in% c("species", "genus"))
  cap <- r$qid[r$rule == "obitag_floor_ident+gap" & !r$rank %in% "species"]
  sp_ref <- r$qid[r$rule == "obitag_floor_ident" & r$rank %in% "species"]
  capped <- intersect(sp_ref, cap)
  flagged_sp <- m$qid[!is.na(m$`Species flag`) & m$`LCA ranked rank` %in% "species"]
  data.frame(
    scenario = sc, queries = nrow(m),
    same_call = round(100 * mean(same_call), 2),
    same_n_tied = round(100 * mean(m$n_tied == m$`N tied`, na.rm = TRUE), 2),
    capped_bench = length(capped), flagged_package = length(flagged_sp),
    both = length(intersect(capped, flagged_sp))
  )
}
res <- rbind(check(s, "in"), check(s[!(s$species %in% s$true_species), ], "out"))
print(res, row.names = FALSE)
