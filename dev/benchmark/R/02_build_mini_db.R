# Build a small multi-marker BLAST DB from core_nt subjects already hit by the
# source projects (all markers), for phase 0 (fast, results comparable).
# Run from repo root: Rscript dev/benchmark/R/02_build_mini_db.R
suppressPackageStartupMessages(library(DBI))
set.seed(20261119)
prj <- path.expand("~/prjcts")
bin <- path.expand("~/.local/share/R/condathis/envs/blastr-blast-env/bin")
core_nt <- "/data/databases/core_nt_nov26/core_nt"
db_dir <- "dev/benchmark/db"
dir.create(db_dir, showWarnings = FALSE)

fs <- c(
  Sys.glob(file.path(prj, "mtbrcdR/analyses/*/db/*.sqlite")),
  Sys.glob(file.path(prj, "mtbrcdR/analyses/*/*/db/*.sqlite")),
  Sys.glob(file.path(prj, "curatR/db/*.data.sqlite"))
)
acc <- lapply(fs, function(f) {
  con <- dbConnect(RSQLite::SQLite(), f, flags = RSQLite::SQLITE_RO)
  on.exit(dbDisconnect(con))
  x <- unique(dbGetQuery(con, "select distinct subject from hits")$subject)
  x <- sub("^gi\\|[0-9]+\\|[a-z]+\\|([^|]+)\\|.*$", "\\1", x)
  x <- x[!grepl("^LGC_", x) & !is.na(x)]
  if (length(x) > 20000) x <- sample(x, 20000)
  x
})
acc <- unique(unlist(acc))
cat("Accessions:", length(acc), "\n")
acc_file <- file.path(db_dir, "mini_acc.txt")
writeLines(acc, acc_file)

# lengths + taxids; keep entries <= 50 kb so whole genomes do not dominate
meta <- system2(file.path(bin, "blastdbcmd"),
  c("-db", core_nt, "-entry_batch", acc_file, "-outfmt", shQuote("%a\t%l\t%T")),
  stdout = TRUE, stderr = FALSE
)
meta <- utils::read.delim(text = meta, header = FALSE, col.names = c("acc", "len", "taxid"))
meta <- meta[!duplicated(meta$acc) & meta$len <= 50000, ]
cat("Kept:", nrow(meta), "entries,", round(sum(meta$len) / 1e6, 1), "Mb\n")
writeLines(meta$acc, acc_file)

fasta <- file.path(db_dir, "mini_multimarker.fasta")
system2(file.path(bin, "blastdbcmd"),
  c("-db", core_nt, "-entry_batch", acc_file, "-outfmt", "%f", "-out", fasta)
)
# blastdbcmd may return an accession more than once: keep the first record
fa <- readLines(fasta)
h <- grepl("^>", fa)
rec <- cumsum(h)
acc_h <- sub("^>([^ ]+).*", "\\1", fa[h])
keep_rec <- which(!duplicated(acc_h))
writeLines(fa[rec %in% keep_rec], fasta)
cat("Unique FASTA records:", length(keep_rec), "\n")
# headers become ">ACC title"; taxid map uses the accession
taxmap <- file.path(db_dir, "mini_multimarker.taxid_map.txt")
utils::write.table(meta[, c("acc", "taxid")], taxmap, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
system2(file.path(bin, "makeblastdb"), c(
  "-in", fasta, "-dbtype", "nucl", "-parse_seqids", "-blastdb_version", "5",
  "-taxid_map", taxmap, "-out", file.path(db_dir, "mini_multimarker"),
  "-title", "mini_multimarker"
))
