# Calibration of the hit ranking (rank_by): which criterion picks the right
# taxon? Public test on core_nt:
#   1. the amplicon sequences of the benchmark pool (12S, 16S, COI, ITS2,
#      rbcL, trnL) are used only to locate the marker regions in core_nt;
#   2. for each marker, up to 300 references of distinct species (at most 3
#      per genus) are cut to the region found (qcov >= 90) -> queries of
#      known taxonomy (public GenBank sequences);
#   3. the queries are searched in core_nt without their own accessions
#      (-negative_seqidlist), qcov >= 50, 100 subjects;
#   4. the hits are ranked offline by identity, by (w * identity + qcov) /
#      (w + 1) for several w (w = 10 is the Similarity score) and by
#      bitscore; the first hit is compared with the true species / genus.
# Results (aggregated) in results/rank_calibration_*.tsv; nothing of the
# client data is written.
# Run from repo root:
#   nohup nice -n 5 Rscript dev/benchmark/R/10_rank_calibration.R > dev/benchmark/results/rank_calibration.log 2>&1 &
suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(10)

db <- "/data/databases/core_nt_nov26/core_nt"
cores <- 32L
out_dir <- "dev/benchmark/results/rank_calibration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
bin <- dirname(blast_bin_path("blastn", "blastr-blast-env"))
markers <- c("12S", "16S", "COI", "ITS2", "rbcL", "trnL")
n_per_marker <- 300L
max_per_genus <- 3L

# 1. locate the marker regions --------------------------------------------
pool <- read.delim("dev/benchmark/data/seqs_pool.tsv")
fa <- readLines("dev/benchmark/data/seqs_pool.fasta")
pool$sequence <- fa[match(paste0(">", pool$seq_id), fa) + 1L]
pool <- pool[pool$marker %in% markers & !is.na(pool$sequence), ]
regions_rds <- file.path(out_dir, "regions.rds")
if (!file.exists(regions_rds)) {
  h <- parallel_blast_hsps(
    pool$sequence, db,
    total_cores = cores, max_target_seqs = 500L, max_hsps = 1L, perc_id = 75, perc_qcov_hsp = 90,
    fields = c("saccver", "pident", "qcovhsp", "sstart", "send", "sstrand", "staxid"),
    verbose = "silent"
  )
  h$marker <- pool$marker[match(h$Sequence, pool$sequence)]
  saveRDS(h, regions_rds)
}
h <- readRDS(regions_rds)
h <- h[!duplicated(h$saccver), ]
message(format(Sys.time()), " regions: ", nrow(h), " subjects")

# 2. references of known species -------------------------------------------
tax <- suppressWarnings(get_tax(h$staxid))
h <- dplyr::left_join(h, tax[, c("taxid", "rank", "species", "genus", "family")], by = c(staxid = "taxid"))
bad <- "( sp\\.| cf\\.| aff\\.| x |uncultured|environmental|unidentified|hybrid)"
h <- h[h$rank %in% "species" & !is.na(h$genus) & !grepl(bad, h$species), ]
h <- h[sample(nrow(h)), ]
h <- h[!duplicated(paste(h$marker, h$species)), ]
h <- h[stats::ave(seq_len(nrow(h)), paste(h$marker, h$genus), FUN = seq_along) <= max_per_genus, ]
ref <- do.call(rbind, lapply(split(h, h$marker), function(x) utils::head(x, n_per_marker)))
message(format(Sys.time()), " references: ", paste(names(table(ref$marker)), table(ref$marker), collapse = ", "))

# cut each reference to its region
eb <- file.path(out_dir, "entry_batch.txt")
writeLines(sprintf(
  "%s %d-%d %s", ref$saccver, pmin(ref$sstart, ref$send), pmax(ref$sstart, ref$send),
  ifelse(ref$sstrand == "minus", "minus", "plus")
), eb)
q <- system2(file.path(bin, "blastdbcmd"), c("-db", db, "-entry_batch", eb, "-outfmt", "%s"), stdout = TRUE)
stopifnot(length(q) == nrow(ref))
ref$query <- q
# identical regions of different species: no single right answer
dup <- ref$query[duplicated(ref$query)]
ref <- ref[!ref$query %in% dup, ]

# 3. search without the references themselves -------------------------------
neg_txt <- file.path(out_dir, "negative.txt")
writeLines(unique(h$saccver[h$saccver %in% ref$saccver]), neg_txt)
neg <- file.path(out_dir, "negative.bsl")
st <- system2(file.path(bin, "blastdb_aliastool"), c("-seqid_file_in", neg_txt, "-seqid_file_out", neg))
stopifnot(st == 0L, file.exists(neg))
hits_rds <- file.path(out_dir, "hits.rds")
if (!file.exists(hits_rds)) {
  s <- parallel_blast_hsps(
    ref$query, db,
    total_cores = cores, max_target_seqs = 100L, max_hsps = 1L, perc_id = 70, perc_qcov_hsp = 50,
    fields = c("saccver", "pident", "length", "qcovhsp", "bitscore", "staxid"),
    extra_args = c("-negative_seqidlist", neg),
    verbose = "silent"
  )
  saveRDS(s, hits_rds)
}
s <- readRDS(hits_rds)
stopifnot(!any(s$saccver %in% ref$saccver))
message(format(Sys.time()), " hits: ", nrow(s), " for ", length(unique(s$Sequence)), " queries")
saveRDS(ref[, c("marker", "saccver", "staxid", "species", "genus", "family", "query")], file.path(out_dir, "references.rds"))
message("Rank calibration searches done")
