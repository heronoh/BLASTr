# Build the stratified seqs pool and nested benchmark subsets.
#
# Sources (read-only): mtbrcdR analysis bases (EM291, EM293, EM294) and curatR
# bases (bee_iDNA, eDNA_Lagoinha, eDNA_Lycengraulis, iDNA_tadpoles).
# Output: dev/benchmark/data/{seqs_pool,S100,S500,S2000}.{fasta,tsv} + MANIFEST.tsv
#
# Run from the BLASTr repo root: Rscript dev/benchmark/R/01_build_seqs.R

suppressPackageStartupMessages({
  library(DBI)
  library(dplyr)
})
set.seed(20261119)

out_dir <- "dev/benchmark/data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prj <- path.expand("~/prjcts")

ro <- function(f) dbConnect(RSQLite::SQLite(), f, flags = RSQLite::SQLITE_RO)

# mtbrcdR bases ---------------------------------------------------------------
read_mtbrcdr <- function(f, source) {
  con <- ro(f)
  on.exit(dbDisconnect(con))
  seqs <- dbGetQuery(con, "select sequence from asvs")
  prim <- dbGetQuery(con, "select distinct sequence, Primer as primer, read_origin from abundance")
  mk <- dbGetQuery(con, "select primer, marker from primers")
  ids <- dbGetQuery(con, "select sequence, pseudo_score from identification")
  seqs |>
    left_join(prim, by = "sequence", relationship = "one-to-many") |>
    left_join(mk, by = "primer") |>
    mutate(marker = ifelse(is.na(marker), primer_marker(primer), marker)) |>
    left_join(ids, by = "sequence") |>
    mutate(source = source)
}

# curatR bases (legacy schema) --------------------------------------------------
read_curatr <- function(f, source) {
  con <- ro(f)
  on.exit(dbDisconnect(con))
  # some bases keep the sequence only in "Sequence (ASV tip)"
  a <- dbGetQuery(con, r"[select "ASV header" as h,
    case when "ASV (Sequence)" in ('-', '') or "ASV (Sequence)" is null
         then "Sequence (ASV tip)" else "ASV (Sequence)" end as sequence,
    "BLASTn pseudo-score" as pseudo_score,
    "Read origin" as read_origin from asvs]")
  if ("Primer" %in% dbListFields(con, "abundance")) {
    p <- dbGetQuery(con, 'select distinct "ASV header" as h, Primer as primer from abundance')
  } else {
    p <- dbGetQuery(con, 'select distinct a."ASV header" as h, l.Primer as primer
                          from abundance a join libraries l using (Unique_File_name)')
  }
  a |>
    left_join(p, by = "h", relationship = "one-to-many") |>
    select(-h) |>
    mutate(marker = primer_marker(primer), source = source)
}

# primer name -> marker, for bases that only store the primer
primer_marker <- function(primer) {
  key <- tolower(primer)
  dplyr::case_when(
    grepl("rbcl", key) ~ "rbcL",
    grepl("trnl", key) ~ "trnL",
    grepl("its", key) ~ "ITS2",
    grepl("mifish|12s|batra", key) ~ "12S",
    grepl("16s", key) ~ "16S",
    grepl("euk|18s", key) ~ "18S",
    grepl("fwh|coi|fish", key) ~ "COI",
    TRUE ~ primer
  )
}

raw <- bind_rows(
  read_mtbrcdr(file.path(prj, "mtbrcdR/analyses/eDNA_vale_dmx/EM291_os1e2/db/dmx_Vale_nov26--complete_results_mtbrcdR.sqlite"), "EM291"),
  read_mtbrcdr(file.path(prj, "mtbrcdR/analyses/MiSeqi100_Seq11/db/EM293--complete_results_mtbrcdR.sqlite"), "EM293"),
  read_mtbrcdr(file.path(prj, "mtbrcdR/analyses/MiSeqi100_Seq11/db/EM294--complete_results_mtbrcdR.sqlite"), "EM294"),
  read_curatr(file.path(prj, "curatR/db/bee_iDNA.data.sqlite"), "bee_iDNA"),
  read_curatr(file.path(prj, "curatR/db/eDNA_Lagoinha.data.sqlite"), "eDNA_Lagoinha"),
  read_curatr(file.path(prj, "curatR/db/eDNA_Lycengraulis.data.sqlite"), "eDNA_Lycengraulis"),
  read_curatr(file.path(prj, "curatR/db/iDNA_tadpoles.data.sqlite"), "iDNA_tadpoles")
)

seqs <- raw |>
  filter(!is.na(sequence), nzchar(sequence)) |>
  mutate(
    sequence = toupper(gsub("\\s", "", sequence)),
    pseudo_score = suppressWarnings(as.numeric(pseudo_score)),
    is_concat = grepl("NNNNN", sequence),
    read_origin = ifelse(is_concat, "concat", "merged"),
    length = nchar(gsub("N", "", sequence))
  ) |>
  # one row per sequence: first source/primer wins, best known score
  group_by(sequence) |>
  summarise(
    source = first(source),
    primer = first(primer),
    marker = first(marker),
    read_origin = first(read_origin),
    length = first(length),
    pseudo_score = suppressWarnings(max(pseudo_score, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  mutate(
    pseudo_score = ifelse(is.finite(pseudo_score), pseudo_score, NA_real_),
    hit_class = case_when(
      is.na(pseudo_score) ~ "no_hit",
      pseudo_score >= 98 ~ "ge98",
      TRUE ~ "lt98"
    ),
    length_bin = cut(length, c(0, 100, 200, 400, Inf),
      labels = c("le100", "101-200", "201-400", "gt400")
    ),
    marker = ifelse(is.na(marker), "unknown", marker)
  ) |>
  filter(length >= 30)

cat("Seqs available:", nrow(seqs), "\n")
print(count(seqs, marker, sort = TRUE), n = 50)

# Stratified pool: up to `k` seqs per marker x length_bin x hit_class x read_origin
k <- 12L
pool <- seqs |>
  group_by(marker, length_bin, hit_class, read_origin) |>
  slice_sample(n = k) |>
  ungroup() |>
  slice_sample(prop = 1) |>
  mutate(seq_id = sprintf("q%04d", row_number()))

cat("Pool:", nrow(pool), "\n")

# Nested subsets: S100 is the most balanced (round-robin over strata).
strata_order <- pool |>
  group_by(marker, length_bin, hit_class, read_origin) |>
  mutate(rank_in_stratum = row_number()) |>
  ungroup() |>
  arrange(rank_in_stratum, sample(n())) |>
  mutate(order = row_number())

subset_n <- c(S100 = 100L, S500 = 500L, S2000 = 2000L)
write_set <- function(df, name) {
  df <- arrange(df, seq_id)
  write.table(
    select(df, seq_id, source, primer, marker, read_origin, length, length_bin, hit_class, pseudo_score),
    file.path(out_dir, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE
  )
  writeLines(paste0(">", df$seq_id, "\n", df$sequence), file.path(out_dir, paste0(name, ".fasta")))
}
write_set(pool, "seqs_pool")
for (nm in names(subset_n)) {
  n <- min(subset_n[[nm]], nrow(pool))
  write_set(filter(strata_order, order <= n), nm)
}

files <- list.files(out_dir, pattern = "\\.(fasta|tsv)$", full.names = TRUE)
write.table(
  data.frame(file = basename(files), md5 = unname(tools::md5sum(files))),
  file.path(out_dir, "MANIFEST.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)

s100 <- read.delim(file.path(out_dir, "S100.tsv"))
cat("\nS100 composition\n")
print(table(s100$marker, s100$hit_class))
print(table(s100$length_bin, s100$read_origin))
