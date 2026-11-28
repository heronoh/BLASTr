# Multi-marker benchmark: resolution of each primer (review of hub
# 2026-11-24-2335, item 3). For each query, with the hits as parallel_blast()
# keeps them by default (identity >= 80, 100 subjects, informative references
# with a named species), whether the references identical to the query, or
# within 0.5% of it, belong to two or more species. Distance to the query:
# differences in the alignment plus query bases outside it, over the query
# length (as the tie rules). Template entry left out (scenario "in").
# Output: results/primer_bench_resolution.tsv, one row per primer and group,
# with the difference in species wrong (each rule minus the current one) out of
# the database, without cut and with the cut of 98.
# Run from repo root: Rscript dev/benchmark/R/22_primer_resolution.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(dplyr))
source("dev/benchmark/R/tie_rules_lib.R")

out_dir <- "dev/benchmark/results/primer_bench"
res_dir <- "dev/benchmark/results"
tax_db <- path.expand("~/.cache/R/BLASTr/taxonomy/taxonomy_2026-11-12.sqlite") # taxonomy of core_nt_nov26
b <- bench_truth(readRDS(file.path(out_dir, "queries.rds")))

s <- readRDS(file.path(out_dir, "hits.rds"))
s <- s[s$Sequence %in% b$truth$seq, c("Sequence", "saccver", "stitle", "pident", "length", "nident", "qstart", "qend", "qlen", "staxid")]
num <- c("pident", "length", "nident", "qstart", "qend", "qlen")
s[num] <- lapply(s[num], as.numeric)
s$qid <- match(s$Sequence, b$truth$seq)
s <- s[s$pident >= 80, ]
s <- s[stats::ave(seq_len(nrow(s)), s$qid, FUN = seq_along) <= 100L, ] # BLAST order, as default_setting()
s$staxid <- trimws(sub(";.*", "", as.character(s$staxid)))
tax <- suppressWarnings(get_tax(unique(s$staxid), taxonomy_db = tax_db))
s$species <- tax$species[match(s$staxid, tax$taxid)]
s <- s[!is.na(s$species) & !grepl(open_nomenclature_pattern, s$species, perl = TRUE), ]
s <- s[is.na(uninformative_category(s$saccver, s$stitle, uninformative_references())), ]
s$dq <- ((s$length - s$nident) + (s$qlen - (s$qend - s$qstart + 1))) / s$qlen

per_query <- s |>
  group_by(qid) |>
  summarise(
    sp_ident = n_distinct(species[dq == 0]),
    sp_05 = n_distinct(species[dq <= 0.005]),
    .groups = "drop"
  )
per_query$seq <- b$truth$seq[per_query$qid]
x <- inner_join(b$qq, per_query, by = "seq")
x$len <- nchar(x$seq)
res <- x |>
  group_by(primer, marker, group) |>
  summarise(
    n = n(), len = median(len),
    ident_any = round(100 * mean(sp_ident >= 1), 1),
    ident_2sp = round(100 * mean(sp_ident >= 2), 1),
    within05_2sp = round(100 * mean(sp_05 >= 2), 1),
    .groups = "drop"
  )

# difference in species wrong out of the database, each rule minus the current one
rules <- c("obitag_floor", "obitag_floor_ident", "obitag_floor_ident+gap")
for (cn in c("", "_cut98")) {
  o <- read.delim(file.path(res_dir, paste0("primer_bench_summary", cn, ".tsv")))
  o <- o[o$scenario == "out", ]
  k <- paste(res$primer, res$group)
  cur <- o$species_wrong[o$rule == "current"][match(k, o$set[o$rule == "current"])]
  for (r in rules) {
    v <- o$species_wrong[o$rule == r][match(k, o$set[o$rule == r])]
    res[[paste0("d_out_", r, cn)]] <- round(v - cur, 1)
  }
}
res <- res[order(-res$within05_2sp), ]
write.table(res, file.path(res_dir, "primer_bench_resolution.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ok <- res$n >= 50
for (v in grep("^d_out_", names(res), value = TRUE)) {
  message(v, ": Spearman with within05_2sp = ", round(stats::cor(res$within05_2sp[ok], res[[v]][ok], method = "spearman", use = "complete"), 2),
          "; with length = ", round(stats::cor(res$len[ok], res[[v]][ok], method = "spearman", use = "complete"), 2))
}
print(as.data.frame(res[ok, ]), row.names = FALSE)
