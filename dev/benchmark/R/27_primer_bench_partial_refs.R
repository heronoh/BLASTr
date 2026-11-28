# Multi-marker benchmark with references that cover only part of the query
# (review of hub 2026-11-25-0218, in 2026-11-25-0222): in the benchmark the
# best hit covers the whole query in 92% of the queries (82% with the species
# out), so the two windows of the rule by references (whole distance, or on
# the aligned region) barely differ there. Added after the real runs.
#
# Each hit, with probability 1/2 (seed 27), loses the query bases beyond a
# coverage `c` of the query, from one end chosen at random (as references
# shorter than the amplicon); its differences, identity, length, nident and
# qcovhsp are recomputed on what is left (identity >= 80 and the first 100
# subjects are then taken as parallel_blast() does). c = 1 is the benchmark as
# it is.
#
# Rules: current and the rule by references (obitag_floor_ident) with the
# weight of the query bases a hit does not cover varied in the choice of the
# best hits and in the window (review, hub 2026-11-25-0244): refs_b<b>_w<w>,
# b and w in 1 (the whole distance, as delivered: refs_b1_w1), 11 (1/11) and 0
# (the rate on the aligned region); b = S chooses the best hits by the
# Similarity score (the rank of BLASTr). Calls cut as mtbrcdR reports them
# (cut98), the cut by the Similarity score of the best informative hit.
# Output (results/):
#   primer_bench_partial_summary.tsv   % species right / wrong, genus right /
#                                      wrong by coverage c, scenario and rule
#   primer_bench_partial_by_cov.tsv    the same by coverage of the best hit
#                                      (the closest by the whole distance)
#   primer_bench_partial_by_marker.tsv the same by marker and group
# Per-query results: results/primer_bench/rules_res_partial_<c>_<group>.rds
# (8 processes: 4 coverages x 2 groups of rules)
# Run from repo root: Rscript dev/benchmark/R/27_primer_bench_partial_refs.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(dplyr))
source("dev/benchmark/R/tie_rules_lib.R")

out_dir <- "dev/benchmark/results/primer_bench"
res_dir <- "dev/benchmark/results"
tax_db <- path.expand("~/.cache/R/BLASTr/taxonomy/taxonomy_2026-11-12.sqlite") # taxonomy of core_nt_nov26
rules <- c("current", as.vector(outer(c("refs_b1", "refs_b11", "refs_bS", "refs_b0"), c("_w1", "_w11", "_w0"), paste0)))
groups <- list(rules[1:7], rules[8:13])
coverages <- c(1, 0.97, 0.9, 0.8)
cut98 <- c(species = 98, genus = 95)

b <- bench_truth(readRDS(file.path(out_dir, "queries.rds")))
# TEST_N: a quick check on the first queries, written under TEST_DIR
test_n <- as.integer(Sys.getenv("TEST_N", "0"))
if (test_n > 0) {
  b$truth <- b$truth[seq_len(test_n), ]
  b$qq <- b$qq[b$qq$seq %in% b$truth$seq, ]
  out_dir <- res_dir <- Sys.getenv("TEST_DIR")
  stopifnot(dir.exists(out_dir))
  coverages <- c(1, 0.9)
}
truth <- b$truth
raw <- readRDS(file.path("dev/benchmark/results/primer_bench", "hits.rds"))
raw <- raw[raw$Sequence %in% truth$seq, ]
num <- c("pident", "length", "nident", "qstart", "qend", "qlen", "qcovhsp", "bitscore")
raw[num] <- lapply(raw[num], as.numeric)
raw$staxid <- trimws(sub(";.*", "", as.character(raw$staxid)))
tax <- suppressWarnings(get_tax(unique(raw$staxid), taxonomy_db = tax_db))

# Each hit, with probability 1/2, keeps only the query region lo..hi (one end
# cut at a coverage c of the query); its differences (`states`), nident,
# length, identity, qcovhsp and distances are recomputed on what is left.
truncate_hits <- function(s, c) {
  if (c >= 1) return(s)
  set.seed(27)
  cut <- stats::runif(nrow(s)) < 0.5
  right <- stats::runif(nrow(s)) < 0.5
  keep_len <- floor(c * s$qlen)
  lo <- ifelse(cut & !right, s$qlen - keep_len + 1, 1)
  hi <- ifelse(cut & right, keep_len, s$qlen)
  nq0 <- pmax(s$qstart, lo)
  nq1 <- pmin(s$qend, hi)
  ok <- nq1 >= nq0
  s <- s[ok, ]
  nq0 <- nq0[ok]
  nq1 <- nq1[ok]
  i <- which(nq0 > s$qstart | nq1 < s$qend)
  if (length(i)) {
    st <- mapply(function(v, a, z) {
      if (!length(v)) return(v)
      nm <- names(v)
      p <- as.integer(substring(nm, 2L))
      v[ifelse(startsWith(nm, "i"), p > a & p <= z, p >= a & p <= z)]
    }, s$states[i], nq0[i], nq1[i], SIMPLIFY = FALSE, USE.NAMES = FALSE)
    n_sub <- vapply(st, function(v) sum(startsWith(names(v), "p")), numeric(1L))
    n_ins <- vapply(st, function(v) sum(nchar(v[startsWith(names(v), "i")])), numeric(1L))
    span <- nq1[i] - nq0[i] + 1
    s$states[i] <- st
    s$nident[i] <- span - n_sub
    s$length[i] <- span + n_ins
    s$qstart[i] <- nq0[i]
    s$qend[i] <- nq1[i]
    s$pident[i] <- 100 * s$nident[i] / s$length[i]
    s$qcovhsp[i] <- round(100 * span / s$qlen[i])
    s$d_aln[i] <- s$length[i] - s$nident[i]
    s$d_cov[i] <- s$d_aln[i] + (s$qlen[i] - span)
    s$diffs[i] <- list(NULL) # not used by these rules
  }
  s
}

cov_band <- function(v) {
  cut(v, c(-Inf, 0.8, 0.9, 0.97, 0.99999, Inf), right = FALSE,
      labels = c("<80", "80-90", "90-97", "97-<100", "100"))
}
# coverage of the best hit (closest by the whole distance, among the ranked
# ones) in each scenario
best_cov <- function(s) {
  bc <- function(h) {
    h <- default_setting(h)
    h <- h[h$informative, ]
    h$dq <- h$d_cov / h$qlen
    h <- h[order(h$qid, h$dq), ]
    h <- h[!duplicated(h$qid), ]
    data.frame(qid = h$qid, best_cov = (h$qend - h$qstart + 1) / h$qlen)
  }
  list(`in` = bc(s), out = bc(s[!(s$species %in% s$true_species), ]))
}

prepped <- prep_hits(raw, match(raw$Sequence, truth$seq), tax, truth$species)
rm(raw)
jobs <- expand.grid(cv = coverages, g = seq_along(groups))
one <- function(k) {
  cv <- jobs$cv[k]
  g <- jobs$g[k]
  f <- file.path(out_dir, sprintf("rules_res_partial_%s_%d.rds", cv, g))
  if (file.exists(f)) return(readRDS(f))
  s <- truncate_hits(prepped, cv)
  twin <- unique(s$qid[(s$species == s$true_species) %in% TRUE & s$nident == s$qlen & s$length == s$qlen])
  r <- run_scenarios(s, tax, tax_db, default_setting, 100L, groups[[g]])
  bc <- best_cov(s)
  runs <- r$runs
  runs$in_noident <- runs[["in"]][!runs[["in"]]$qid %in% twin, ]
  runs$in_noident$scenario <- "in_noident"
  for (sc in names(runs)) {
    bsc <- if (sc == "out") bc$out else bc[["in"]]
    runs[[sc]] <- merge(runs[[sc]], bsc, by = "qid", all.x = TRUE)
  }
  x <- do.call(rbind, runs)
  x$coverage <- cv
  saveRDS(x, f)
  message(format(Sys.time()), " coverage ", cv, " group ", g, " done")
  x
}
all <- parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = nrow(jobs), mc.preschedule = FALSE)
bad <- vapply(all, inherits, logical(1L), "try-error")
if (any(bad)) stop("failed jobs: ", paste(which(bad), collapse = ", "), "\n", all[bad][[1L]])
rm(prepped)

# cut98 as mtbrcdR reports the calls: by the Similarity score of the best hit
# tied by each rule (mtbrcdR 657887b, "tied"), or of the best hit of all, as
# the earlier runs of the benchmark (19_, "hit1")
x0 <- do.call(rbind, all)
cut_calls <- function(x, score) {
  x$true_species <- truth$species[x$qid]
  x$true_genus <- truth$genus[x$qid]
  k <- x$rank %in% "species" & !(score >= cut98[["species"]]) %in% TRUE
  x$rank[k] <- "genus"
  x$species[k] <- NA_character_
  k <- x$rank %in% c("species", "genus") & !(score >= cut98[["genus"]]) %in% TRUE
  x$rank[k] <- "above genus"
  x$species[k] <- NA_character_
  x$genus[k] <- NA_character_
  k <- x$rank %in% "species" & grepl(open_nomenclature_pattern, x$species, perl = TRUE)
  x$rank[k] <- "genus"
  x$species[k] <- NA_character_
  x$sp_ok <- x$rank %in% "species" & (x$species == x$true_species) %in% TRUE
  x$sp_wrong <- x$rank %in% "species" & !x$sp_ok
  x$gen_ok <- x$rank %in% c("species", "genus") & (x$genus == x$true_genus) %in% TRUE
  x$gen_wrong <- x$rank %in% c("species", "genus") & !x$gen_ok
  x$cov_band <- cov_band(x$best_cov)
  x
}
x <- rbind(cbind(cut_score = "tied", cut_calls(x0, x0$score_tied)),
           cbind(cut_score = "hit1", cut_calls(x0, x0$score)))
pct <- function(v) round(100 * mean(v), 1)
summ <- function(d, by) {
  d |>
    group_by(across(all_of(by))) |>
    summarise(n = n(), species_right = pct(sp_ok), species_wrong = pct(sp_wrong),
              genus_right = pct(gen_ok), genus_wrong = pct(gen_wrong), median_tied = median(n_tied),
              .groups = "drop") |>
    arrange(across(all_of(setdiff(by, "rule"))), match(rule, rule_order))
}
# one row per primer and query, as in 19_
x$seq <- truth$seq[x$qid]
x <- inner_join(x, b$qq, by = "seq", relationship = "many-to-many")
x$marker_group <- paste(x$marker, x$group)
s1 <- summ(x, c("cut_score", "coverage", "scenario", "rule"))
x <- x[x$cut_score == "tied", ]
s2 <- summ(x, c("coverage", "scenario", "cov_band", "rule"))
s3 <- summ(x, c("coverage", "scenario", "marker_group", "rule"))
write.table(s1, file.path(res_dir, "primer_bench_partial_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(s2, file.path(res_dir, "primer_bench_partial_by_cov.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(s3, file.path(res_dir, "primer_bench_partial_by_marker.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(s1), row.names = FALSE)
message("finished")
