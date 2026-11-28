# A reference point for the fraction f of sequences whose species is absent
# from the reference database (decision procedure of the tie rules, step 2):
# in real analyses, the fraction of sequences whose best hit is under 97, 98 or
# 99% identity, among those with a best hit of at least 80% identity and 90%
# coverage, by group (vertebrates, invertebrates, plants) from the stored
# taxonomy. Under the cut it mixes species absent from the database with
# variation within the species and sequencing errors: an estimate, not a
# measure.
# Sequences (client data, read only): curatR data bases in BLASTR_BENCH_SEQ_DBS and
# parallel_blast() results saved as .rds in BLASTR_BENCH_BLAST_RES (both
# ":"-separated); only aggregated numbers are written
# (results/absent_fraction.tsv), without the names of the projects.
# Run from repo root:
#   BLASTR_BENCH_SEQ_DBS=<a.sqlite>:<b.sqlite> BLASTR_BENCH_BLAST_RES=<c.rds> \
#     Rscript dev/benchmark/R/20_absent_fraction.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(dplyr))
split_env <- function(v) Filter(nzchar, strsplit(Sys.getenv(v), ":", fixed = TRUE)[[1]])
dbs <- split_env("BLASTR_BENCH_SEQ_DBS")
rds <- split_env("BLASTR_BENCH_BLAST_RES")
if (!length(dbs) && !length(rds)) stop("Set BLASTR_BENCH_SEQ_DBS or BLASTR_BENCH_BLAST_RES.")
tax_db <- path.expand("~/.cache/R/BLASTr/taxonomy/taxonomy_2026-11-12.sqlite")
to_group <- function(kingdom, phylum) {
  ifelse(kingdom %in% "Metazoa" & phylum %in% "Chordata", "vertebrates",
    ifelse(kingdom %in% "Metazoa", "invertebrates", ifelse(kingdom %in% "Viridiplantae", "plants", NA))
  )
}

one <- function(f, i) {
  con <- DBI::dbConnect(RSQLite::SQLite(), f, flags = RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(con))
  q <- function(sql) DBI::dbGetQuery(con, sql)
  h1 <- q('select "ASV header" as h, indentity as id, qcovhsp as cov from hits where cast(hit as integer) = 1')
  tx <- q('select "ASV header" as h, "Kingdom (NCBI)" as kingdom, "Phylum (NCBI)" as phylum from taxonomy')
  d <- merge(h1, tx, by = "h")
  d$id <- suppressWarnings(as.numeric(d$id))
  d$cov <- suppressWarnings(as.numeric(d$cov))
  d$group <- to_group(d$kingdom, d$phylum)
  d <- d[!is.na(d$group) & !is.na(d$id) & d$id >= 80 & d$cov >= 90, ]
  d$dataset <- paste0("dataset", i)
  d[, c("dataset", "group", "id")]
}
# parallel_blast() output: best identity of each sequence among hits with 90%
# coverage, group from the taxid of that hit
one_rds <- function(f, i) {
  x <- as.data.frame(readRDS(f))
  x <- x[!is.na(x$identity) & x$qcovhsp >= 90, c("Sequence", "identity", "staxid")]
  x <- x[order(x$Sequence, -x$identity), ]
  x <- x[!duplicated(x$Sequence), ]
  x$staxid <- trimws(sub(";.*", "", as.character(x$staxid)))
  tx <- suppressWarnings(get_tax(unique(x$staxid), taxonomy_db = tax_db))
  k <- match(x$staxid, tx$taxid)
  d <- data.frame(group = to_group(tx$kingdom[k], tx$phylum[k]), id = x$identity)
  d <- d[!is.na(d$group) & d$id >= 80, ]
  d$dataset <- paste0("dataset", length(dbs) + i)
  d
}
d <- bind_rows(
  lapply(seq_along(dbs), function(i) one(dbs[i], i)),
  lapply(seq_along(rds), function(i) one_rds(rds[i], i))
)
summ <- function(x) {
  x |>
    summarise(
      seqs = n(),
      below_97 = round(100 * mean(id < 97), 1),
      below_98 = round(100 * mean(id < 98), 1),
      below_99 = round(100 * mean(id < 99), 1),
      .groups = "drop"
    )
}
out <- bind_rows(
  summ(group_by(d, dataset, group)),
  mutate(summ(group_by(d, group)), dataset = "pooled")
)
write.table(out, "dev/benchmark/results/absent_fraction.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(out))
