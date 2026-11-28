# Calibration of the rules that choose the hits tied with the best one
# (`at_best`, the voters of blast_lca()). Proposal of mtbrcdR
# (hub 2026-11-24-1907) and the review's requests (2026-11-24-2010).
#
# Queries: the public references of 10_rank_calibration.R (GenBank sequences
# of known species cut to the marker region), searched by 14_tie_rules_search.R
# in core_nt without their own accessions (500 subjects).
# Two scenarios: the species of the query is in the database ("in": other
# sequences of it remain) or not ("out": all the hits of its species removed).
# Main setting as parallel_blast() by default: identity >= 80, the first 100
# subjects in BLAST order (search_alignments); also with all 500.
#
# Rules for the tied set (always among the informative hits, ranked by the
# Similarity score as today):
#   current       Similarity score >= best - 0.5 (BLASTr today)
#   aln_k0/k1     d = length - nident (columns of the alignment that differ);
#                 tied if d <= min(d) + k
#   cov_k0/k1     d as above + the query bases outside the alignment
#   common_k0     compared with the first hit only on the query region both
#                 cover (btop); tied if d_hit <= d_first + k, k = 0
#   common_prop   same, k = floor(0.005 * length of the common region)
#                 (0 bases up to 199 bp, 3 in a COI of 652)
#   (common region shorter than half the query: tied)
# Variants:
#   +sp       a hit in open nomenclature within a genus ("Leporinus sp. 3")
#             votes only when it is strictly closer than every named species
#             tied with it (in the metric of the rule)
#   +voucher  a hit in open nomenclature with the voucher code of a named hit
#             tied with it, at the same distance, does not vote (counted once,
#             with the formal name)
#
# For each rule: the ancestor of blast_lca() (consensus 0.9, named_first) is
# compared with the true species and genus. Writes only aggregated numbers
# (results/tie_rules_summary.tsv, tie_rules_voucher.tsv); public data only.
# Run from repo root: Rscript dev/benchmark/R/15_tie_rules.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(dplyr))
source("dev/benchmark/R/tie_rules_lib.R")

d_dir <- "dev/benchmark/results/tie_rules"
tax_db <- path.expand("~/.cache/R/BLASTr/taxonomy/taxonomy_2026-11-12.sqlite") # taxonomy of core_nt_nov26
stopifnot(file.exists(tax_db))
ref <- readRDS("dev/benchmark/results/rank_calibration/references.rds")
s <- readRDS(file.path(d_dir, "hits.rds"))
s$staxid <- trimws(sub(";.*", "", as.character(s$staxid)))
tax <- suppressWarnings(get_tax(unique(s$staxid), taxonomy_db = tax_db))
s <- prep_hits(s, match(s$Sequence, ref$query), tax, ref$species)

run <- function(scenario, setting) {
  h <- s
  if (scenario == "out") h <- h[!(h$species %in% h$true_species), ]
  if (setting == "default") h <- default_setting(h)
  res <- run_rules(h, tax, tax_db)
  res$scenario <- scenario
  res$setting <- setting
  res
}

t0 <- Sys.time()
runs <- list()
for (sc in c("in", "out")) for (st in c("default", "500")) {
  runs[[paste(sc, st)]] <- run(sc, st)
  message(format(Sys.time()), " ", sc, " ", st, " done")
}
res <- do.call(rbind, runs)
res$marker <- ref$marker[res$qid]
res$true_species <- ref$species[res$qid]
res$true_genus <- ref$genus[res$qid]
qlen <- tapply(nchar(ref$query), ref$marker, stats::median)
res$group <- ifelse(qlen[res$marker] <= 200, "short", "long")
out <- summarise_rules(res, c("marker", "group"))
write.table(out, "dev/benchmark/results/tie_rules_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# voucher codes: how often they are read, and how often the rule acts
vd <- do.call(rbind, lapply(names(runs), function(n) {
  v <- attr(runs[[n]], "voucher_dropped")
  data.frame(run = n, rule = names(v), open_hits_dropped = as.integer(v))
}))
vstat <- data.frame(
  hits = nrow(s), with_voucher = sum(!is.na(s$voucher)),
  open_hits = sum(s$open), open_with_voucher = sum(s$open & !is.na(s$voucher))
)
write.table(vd, "dev/benchmark/results/tie_rules_voucher.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(vstat)
print(vd)
message("done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
