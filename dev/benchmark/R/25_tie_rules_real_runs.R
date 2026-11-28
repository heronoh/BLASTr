# Check of the comparison of the two tie rules in real runs of mtbrcdR (hub
# 2026-11-25-0003, -0006; tables of mtbrcdR in 2026-11-25-0203, -0208), on the
# BLASTr side:
#   1. `at_best_refs` stored in the run, recomputed from its btop;
#   2. the base of mtbrcdR against blast_lca() of the package (both rules,
#      species flag);
#   3. hit 1 in the candidates of the rule by references (mtbrcdR adds it);
#   4. added after seeing the runs: share of the window that comes from query
#      bases the best reference does not cover, and the LCA with the window
#      measured on the aligned region only (best still chosen by the full
#      distance).
# Reads only the run folder (RUN_DIR) and writes aggregated numbers to OUT.
# Run from repo root:
#   RUN_DIR=<run folder> OUT=<dir> Rscript dev/benchmark/R/25_tie_rules_real_runs.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(dplyr))

run <- Sys.getenv("RUN_DIR")
out <- Sys.getenv("OUT")
stopifnot(dir.exists(run), dir.exists(out))
tag <- basename(run)
x <- readRDS(file.path(run, "results/blast/blast_res.rds"))
db <- list.files(file.path(run, "final_results"), "complete_results_mtbrcdR$",
                 recursive = TRUE, include.dirs = TRUE, full.names = TRUE)
idt <- readr::read_tsv(file.path(db, "identification.tsv"), show_col_types = FALSE,
                       col_types = readr::cols(.default = "c"), na = c("", "NA"))
ad <- readr::read_tsv(file.path(run, "comparacao_regras_empate/asvs_detalhe.tsv"), show_col_types = FALSE)
tx <- suppressWarnings(get_tax(x, ranks = c("phylum", "class", "order", "family", "genus", "species")))
x$qlen <- nchar(x$Sequence)
pool <- !is.na(x$rank)
lv <- c(species = 6, genus = 5, family = 4, order = 3, class = 2, phylum = 1)
lca <- function(y, ties = "at_best") {
  suppressWarnings(blast_lca(y, ties = ties, taxonomy = tx))
}
res <- list()

# 1. at_best_refs from btop
rt <- reference_ties(x$Sequence, pool, x$btop, x$`query start`, x$`query end`, x$length, x$nident, x$qlen)
res$recompute <- tibble::tibble(ranked_rows = sum(pool),
                                differences = sum((rt$tied %in% TRUE) != (x$at_best_refs %in% TRUE)))

# 2. base against the package
l0 <- lca(x)
l1 <- lca(x, "at_best_refs")
cmp <- idt |>
  inner_join(l0 |> select(sequence = Sequence, n0 = `N tied`), by = "sequence") |>
  inner_join(l1 |> select(sequence = Sequence, n1 = `N tied`, f1 = `Species flag`, fs1 = `Flagged species`),
             by = "sequence")
eq <- function(a, b) sum(a == b | (is.na(a) & is.na(b)), na.rm = TRUE)
res$base <- tibble::tibble(
  seqs = nrow(cmp), n_tied = eq(cmp$n_tied, as.character(cmp$n0)), n_tied_refs = eq(cmp$n_tied_refs, as.character(cmp$n1)),
  species_flag = eq(cmp$species_flag, cmp$f1), flagged_species = eq(cmp$flagged_species, cmp$fs1)
)

# 3. hit 1 (first ranked) not tied by references
h1 <- x |> filter(pool) |> group_by(Sequence) |> slice_min(rank, n = 1, with_ties = FALSE) |> ungroup()
off <- h1$Sequence[h1$at_best_refs %in% FALSE]
x2 <- x
x2$at_best_refs[x2$Sequence %in% off & paste(x2$Sequence, x2$subject, x2$rank) %in%
                  paste(h1$Sequence, h1$subject, h1$rank)] <- TRUE
l2 <- lca(x2, "at_best_refs")
j3 <- tibble::tibble(Sequence = off) |>
  left_join(l1 |> select(Sequence, s1 = `LCA ranked name`, k1 = `LCA ranked rank`), by = "Sequence") |>
  left_join(l2 |> select(Sequence, s2 = `LCA ranked name`, k2 = `LCA ranked rank`), by = "Sequence") |>
  left_join(ad |> select(Sequence = sequence, categoria, tamanho, reads), by = "Sequence") |>
  mutate(changed = !(s1 == s2) %in% TRUE)
res$hit1 <- j3 |> count(changed, k1, k2, categoria, tamanho, wt = NULL, name = "seqs") |>
  left_join(j3 |> group_by(changed, k1, k2, categoria, tamanho) |> summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop"),
            by = c("changed", "k1", "k2", "categoria", "tamanho"))

# 4. window on the aligned region (added after seeing the runs)
ties_aln <- function(query, pool, btop, qstart, qend, length, nident, qlen, floor = refs_floor_default,
                     min_common = refs_min_common) {
  n <- base::length(query)
  tied <- rep(FALSE, n)
  d <- query_distance(length, nident, qstart, qend, qlen)
  da <- (length - nident) / (qend - qstart + 1)
  pool <- pool %in% TRUE & !is.na(d)
  for (rows in split(seq_len(n), factor(query, levels = unique(query)))) {
    rows <- rows[pool[rows]]
    if (!base::length(rows)) next
    if (anyNA(btop[rows])) {
      tied[rows] <- NA
      next
    }
    st <- btop_states(btop[rows], qstart[rows])
    b <- d[rows] <= min(d[rows]) + 1e-12
    t <- max(min(da[rows][b]), floor)
    sig <- vapply(st[b], function(v) paste(names(v), v, collapse = ";"), character(1L))
    inc <- b
    for (i in which(b)[!duplicated(sig)]) {
      for (k in which(!inc)) {
        lo <- max(qstart[rows[i]], qstart[rows[k]])
        hi <- min(qend[rows[i]], qend[rows[k]])
        L <- hi - lo + 1
        inc[k] <- L < min_common * qlen[rows[i]] || state_diff(st[[i]], st[[k]], lo, hi) / L <= t + 1e-9
      }
    }
    tied[rows] <- inc
  }
  tied
}
xa <- x
xa$at_best_refs <- ifelse(pool, ties_aln(x$Sequence, pool, x$btop, x$`query start`, x$`query end`,
                                         x$length, x$nident, x$qlen), NA)
la <- lca(xa, "at_best_refs")
w <- x |> filter(pool) |>
  mutate(d = query_distance(length, nident, `query start`, `query end`, qlen),
         unc = 1 - (`query end` - `query start` + 1) / qlen) |>
  group_by(Sequence) |> slice_min(d, n = 1, with_ties = FALSE) |> ungroup() |> select(Sequence, d, unc)
j4 <- ad |> select(Sequence = sequence, categoria, tamanho, faixa, reads) |>
  inner_join(w, by = "Sequence") |>
  left_join(l0 |> select(Sequence, n0 = `LCA ranked name`, r0 = `LCA ranked rank`), by = "Sequence") |>
  left_join(l1 |> select(Sequence, n1 = `LCA ranked name`, r1 = `LCA ranked rank`), by = "Sequence") |>
  left_join(la |> select(Sequence, n2 = `LCA ranked name`, r2 = `LCA ranked rank`), by = "Sequence")
dir <- function(a, b) {
  a <- lv[a]; b <- lv[b]
  ifelse(is.na(a) | is.na(b), NA, ifelse(a > b, "current deeper", ifelse(a < b, "refs deeper", "same rank")))
}
res$window_share <- j4 |> filter(d > refs_floor_default) |> group_by(categoria) |>
  summarise(seqs = n(), median_d = median(d), median_uncovered = median(unc), median_share = median(unc / d),
            .groups = "drop")
res$window_aln <- j4 |> group_by(tamanho, refs = dir(r0, r1), refs_aln = dir(r0, r2)) |>
  summarise(seqs = n(), reads = sum(reads), .groups = "drop")
res$window_aln_pairs <- j4 |> filter(tamanho == "in range", !(n0 == n1 & n1 == n2) %in% TRUE) |>
  group_by(n0, n1, n2) |>
  summarise(seqs = n(), reads = sum(reads), median_d = median(d), median_uncovered = median(unc), .groups = "drop") |>
  arrange(desc(reads)) |> head(15)

for (k in names(res)) {
  cat("\n##", k, "\n")
  print(as.data.frame(res[[k]]), row.names = FALSE, digits = 3)
  readr::write_tsv(res[[k]], file.path(out, paste0(tag, "--", k, ".tsv")))
}
