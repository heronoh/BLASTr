# Calibration of the tie rules (see 15_tie_rules.R) on real sequences of several
# markers and of vertebrates, invertebrates and plants (microorganisms left
# out: their references are scarcer and less reliable).
#
# Sequences (client data): from a curatR data base, BLASTR_BENCH_SEQ_DB (read
# only; nothing of it is written outside dev/benchmark/results, which is not
# versioned; only aggregated numbers in results/tie_rules_real_seqs_summary.tsv).
# Selected: sequences whose stored first hit is >= 98% identity and >= 95%
# coverage in one of the three groups, up to 400 per primer and group.
# Search: as 14_tie_rules_search.R (core_nt, 500 subjects, btop).
#
# Truth: the sequences whose identical hits (100% identity over the whole sequence)
# all name the same species (with a name, not "sp."); the group is that of
# the species. Three scenarios:
#   exact  all the hits (the identical references are there)
#   in     the identical hits removed; only sequences whose species still has
#          other hits (the species is in the database, not this haplotype)
#   out    every hit of the species removed (any species called is wrong)
# Run from repo root:
#   BLASTR_CORE_NT_DB=<core_nt prefix> BLASTR_BENCH_SEQ_DB=<curatR data sqlite> \
#     nohup nice -n 5 Rscript dev/benchmark/R/16_tie_rules_real_seqs.R > dev/benchmark/results/tie_rules_real_seqs.log 2>&1 &
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(dplyr))
source("dev/benchmark/R/tie_rules_lib.R")
set.seed(16)

db <- Sys.getenv("BLASTR_CORE_NT_DB")
seq_db <- Sys.getenv("BLASTR_BENCH_SEQ_DB")
if (!nzchar(db) || !nzchar(seq_db)) stop("Set BLASTR_CORE_NT_DB and BLASTR_BENCH_SEQ_DB.")
tax_db <- path.expand("~/.cache/R/BLASTr/taxonomy/taxonomy_2026-11-12.sqlite") # taxonomy of core_nt_nov26
stopifnot(file.exists(tax_db))
cores <- 30L
max_per_cell <- 400L
out_dir <- "dev/benchmark/results/tie_rules_real_seqs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. sequences ---------------------------------------------------------------
seq_rds <- file.path(out_dir, "seqs.rds")
if (!file.exists(seq_rds)) {
  con <- DBI::dbConnect(RSQLite::SQLite(), seq_db, flags = RSQLite::SQLITE_RO)
  q <- function(sql) DBI::dbGetQuery(con, sql)
  # column and table names as the curatR data base has them
  a <- q('select "ASV header" as h, "Sequence (ASV tip)" as seq from asvs')
  p <- q('select distinct "ASV header" as h, Primer as primer from abundance')
  tx <- q('select "ASV header" as h, "Kingdom (NCBI)" as kingdom, "Phylum (NCBI)" as phylum from taxonomy')
  h1 <- q('select "ASV header" as h, indentity as id, qcovhsp as cov from hits where cast(hit as integer) = 1')
  DBI::dbDisconnect(con)
  a$primer <- p$primer[match(a$h, p$h)]
  a <- cbind(a, tx[match(a$h, tx$h), c("kingdom", "phylum")])
  a$id1 <- suppressWarnings(as.numeric(h1$id[match(a$h, h1$h)]))
  a$cov1 <- suppressWarnings(as.numeric(h1$cov[match(a$h, h1$h)]))
  a$group0 <- ifelse(a$kingdom %in% "Metazoa" & a$phylum %in% "Chordata", "vertebrates",
    ifelse(a$kingdom %in% "Metazoa", "invertebrates", ifelse(a$kingdom %in% "Viridiplantae", "plants", NA))
  )
  a$seq <- toupper(gsub("\\s", "", a$seq))
  a <- a[!is.na(a$group0) & !is.na(a$id1) & a$id1 >= 98 & a$cov1 >= 95 &
    !grepl("N", a$seq) & nchar(a$seq) >= 50 & !is.na(a$primer), ]
  a <- a[!duplicated(a$seq), ]
  a <- a[sample(nrow(a)), ]
  a <- a[stats::ave(seq_len(nrow(a)), paste(a$primer, a$group0), FUN = seq_along) <= max_per_cell, ]
  stopifnot(nrow(a) > 0L)
  saveRDS(a[, c("h", "seq", "primer", "group0")], seq_rds)
}
a <- readRDS(seq_rds)
message(format(Sys.time()), " sequences selected: ", nrow(a))
print(table(a$primer, a$group0))

# 2. search -----------------------------------------------------------------
hits_rds <- file.path(out_dir, "hits.rds")
if (!file.exists(hits_rds)) {
  t0 <- Sys.time()
  s <- parallel_blast_hsps(
    a$seq, db,
    total_cores = cores, max_target_seqs = 500L, max_hsps = 1L, perc_id = 70, perc_qcov_hsp = 50,
    fields = c(
      "saccver", "stitle", "pident", "length", "nident", "qstart", "qend", "qlen",
      "qcovhsp", "bitscore", "staxid", "btop"
    ),
    verbose = "silent"
  )
  saveRDS(s, hits_rds)
  message("search: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
}
s <- readRDS(hits_rds)

# 3. truth ------------------------------------------------------------------
s$staxid <- trimws(sub(";.*", "", as.character(s$staxid)))
tax <- suppressWarnings(get_tax(unique(s$staxid), taxonomy_db = tax_db))
qid <- match(s$Sequence, a$seq)
s <- prep_hits(s, qid, tax, rep(NA_character_, nrow(a)))
s$exact <- s$nident == s$qlen & s$length == s$qlen
ex <- s[s$exact & s$informative & !is.na(s$species), ]
truth <- ex |>
  group_by(qid) |>
  summarise(
    n_sp = n_distinct(species), all_named = all(named),
    species = first(species), genus = first(genus), taxid = first(staxid), .groups = "drop"
  ) |>
  filter(n_sp == 1L, all_named)
lin <- tax$lineage_taxids[match(truth$taxid, tax$taxid)]
has <- function(id) grepl(paste0("(^|; )", id, "(;|$)"), lin)
truth$group <- ifelse(has(7742), "vertebrates", ifelse(has(33208), "invertebrates", ifelse(has(33090), "plants", NA)))
truth <- truth[!is.na(truth$group), ]
truth$primer <- a$primer[truth$qid]
message(format(Sys.time()), " sequences with a true species: ", nrow(truth))
s <- s[s$qid %in% truth$qid, ]
s$true_species <- truth$species[match(s$qid, truth$qid)]

# 4. scenarios --------------------------------------------------------------
scen <- list(
  exact = s,
  `in` = local({
    h <- s[!s$exact, ]
    keep <- unique(h$qid[h$species %in% h$true_species & h$informative])
    h[h$qid %in% keep, ]
  }),
  out = s[!(s$species %in% s$true_species), ]
)
runs <- lapply(names(scen), function(sc) {
  r <- run_rules(default_setting(scen[[sc]]), tax, tax_db)
  r$scenario <- sc
  r$setting <- "default"
  message(format(Sys.time()), " ", sc, " done")
  r
})
res <- do.call(rbind, runs)
k <- match(res$qid, truth$qid)
res$true_species <- truth$species[k]
res$true_genus <- truth$genus[k]
res$group <- truth$group[k]
res$primer <- truth$primer[k]
res$primer_group <- paste(res$primer, res$group)
out <- summarise_rules(res, c("group", "primer", "primer_group"))
write.table(out, "dev/benchmark/results/tie_rules_real_seqs_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(table(truth$primer, truth$group))
print(out[out$set %in% c("all", "vertebrates", "invertebrates", "plants") &
  out$rule %in% c("current", "current+sp", "cov_k0", "common_k0", "common_prop", "common_prop+sp"), ], row.names = FALSE)
