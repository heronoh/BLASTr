# Multi-marker benchmark, step 1: queries of known species for every primer
# pair of a primer list (mtbrcdR's tables), by in-silico PCR on public
# GenBank templates of vertebrates, invertebrates and plants.
#
#   1. primer pairs: the tables in BLASTR_BENCH_PRIMERS (":"-separated paths;
#      columns primer, fwd_seq, rev_seq, marker, genome, amplicon_max);
#      prokaryotic markers and primer cocktails (";") left out; 5-base tags
#      of tagged variants ("rbcl-A", "ITS2-B"...) and of primers named "*tag"
#      removed; identical pairs
#      (in either orientation) kept once.
#   2. templates (accessions from NCBI E-utilities, sequences from the local
#      core_nt): complete mitochondrial genomes of vertebrates and of other
#      animals, complete chloroplast genomes of plants, plant and animal
#      rDNA with 18S, 5.8S and 26S/28S (ITS) and complete 18S; named species
#      only, one per genus, up to `n_templates` per kind.
#   3. in-silico PCR (Biostrings, IUPAC codes, up to 2 mismatches per primer
#      under 20 nt and 3 from 20 nt, both strands): the shortest product of
#      each template, without the primers, 25 bp to `max_insert` (or 1.5 x
#      amplicon_max); only A, C, G, T.
#   4. queries: up to `n_per_cell` per primer and group (vertebrates,
#      invertebrates, plants), from distinct templates.
# Writes results/primer_bench/queries.rds and negative.bsl (every accession of
# the template entries, for the search). Public data only.
# Run from repo root:
#   BLASTR_CORE_NT_DB=<core_nt prefix> BLASTR_BENCH_PRIMERS=<primers.tsv>:<more.tsv> \
#     Rscript dev/benchmark/R/17_primer_bench_queries.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(Biostrings))
set.seed(17)

db <- Sys.getenv("BLASTR_CORE_NT_DB")
primer_files <- strsplit(Sys.getenv("BLASTR_BENCH_PRIMERS"), ":", fixed = TRUE)[[1]]
if (!nzchar(db) || !length(primer_files)) stop("Set BLASTR_CORE_NT_DB and BLASTR_BENCH_PRIMERS.")
tax_db <- path.expand("~/.cache/R/BLASTr/taxonomy/taxonomy_2026-11-12.sqlite") # taxonomy of core_nt_nov26
out_dir <- "dev/benchmark/results/primer_bench"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
bin <- dirname(blast_bin_path("blastn", "blastr-blast-env"))
n_templates <- 600L
n_per_cell <- 200L
max_insert <- 1000L

# 1. primer pairs -----------------------------------------------------------
p <- do.call(rbind, lapply(primer_files, function(f) {
  x <- utils::read.delim(f, colClasses = "character")
  x[, c("primer", "fwd_name", "rev_name", "fwd_seq", "rev_seq", "marker", "genome", "amplicon_max", "tag_length")]
}))
p <- p[!duplicated(p$primer), ]
p$fwd_seq <- toupper(gsub("\\s", "", p$fwd_seq))
p$rev_seq <- toupper(gsub("\\s", "", p$rev_seq))
# inosine pairs with any base
p$fwd_seq <- gsub("I", "N", p$fwd_seq, fixed = TRUE)
p$rev_seq <- gsub("I", "N", p$rev_seq, fixed = TRUE)
# MiFish-U as published (Miya et al. 2015; the table lacks one base in each)
mifish <- p$primer == "MiFish"
p$fwd_seq[mifish] <- "GTCGGTAAAACTCGTGCCAGC"
p$rev_seq[mifish] <- "CATAGTGGGGTATCTAATCCCAGTTTG"
p <- p[nzchar(p$fwd_seq) & nzchar(p$rev_seq) & !grepl(";", p$fwd_seq) & !grepl(";", p$rev_seq), ]
p <- p[!p$genome %in% "prokaryotic", ]
# tagged variants: two or more primers named <base>-<letter>
base_name <- sub("-[A-H]$", "", p$primer)
tagged <- grepl("-[A-H]$", p$primer) & stats::ave(seq_along(base_name), base_name, FUN = length) > 1L
# and primers whose names say they carry a tag ("rbcl-357F *tag"), or with tag_length
tl <- suppressWarnings(as.integer(p$tag_length))
named_tag <- grepl("tag", paste(p$fwd_name, p$rev_name), ignore.case = TRUE)
tlen <- ifelse(!is.na(tl), tl, ifelse(tagged | named_tag, 5L, 0L))
p$fwd_seq <- substring(p$fwd_seq, tlen + 1L)
p$rev_seq <- substring(p$rev_seq, tlen + 1L)
p$primer[tagged] <- sub("-[A-H]$", "", p$primer[tagged])
key <- ifelse(p$fwd_seq < p$rev_seq, paste(p$fwd_seq, p$rev_seq), paste(p$rev_seq, p$fwd_seq))
p <- p[!duplicated(key), ]
# kind of template for each primer
p$marker[!nzchar(p$marker) & grepl("ITS", p$primer)] <- "ITS"
p$genome[!nzchar(p$genome) & grepl("ITS", p$primer)] <- "nuclear"
p$kinds <- ifelse(p$genome == "mitochondrial", "vert_mito;invert_mito",
  ifelse(p$genome == "chloroplast", "plant_cp",
    ifelse(p$genome == "nuclear" & grepl("^ITS", p$marker), "plant_its;metazoa_its",
      ifelse(p$genome == "nuclear" & p$marker == "18S", "metazoa_18s;plant_18s", NA)
    )
  )
)
p <- p[!is.na(p$kinds), ]
p$max_len <- ifelse(is.na(suppressWarnings(as.numeric(p$amplicon_max))), max_insert,
  pmin(max_insert, ceiling(1.5 * as.numeric(p$amplicon_max)))
)
message("primer pairs: ", nrow(p))

# 2. templates --------------------------------------------------------------
eutils <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
esearch_acc <- function(term, retmax = 10000L) {
  u <- paste0(
    eutils, "?db=nuccore&idtype=acc&retmax=", retmax, "&tool=blastr_benchmark&term=",
    utils::URLencode(term, reserved = TRUE)
  )
  x <- paste(readLines(u, warn = FALSE), collapse = "")
  Sys.sleep(0.4)
  regmatches(x, gregexpr("(?<=<Id>)[^<]+(?=</Id>)", x, perl = TRUE))[[1]]
}
kinds <- list(
  vert_mito = list(term = 'Vertebrata[Organism] AND mitochondrion[filter] AND "complete genome"[Title] NOT refseq[filter]', groups = "vertebrates"),
  invert_mito = list(term = 'Metazoa[Organism] NOT Vertebrata[Organism] AND mitochondrion[filter] AND "complete genome"[Title] NOT refseq[filter]', groups = "invertebrates"),
  plant_cp = list(term = 'Viridiplantae[Organism] AND chloroplast[filter] AND "complete genome"[Title] NOT refseq[filter]', groups = "plants"),
  plant_its = list(term = "Viridiplantae[Organism] AND 18S[Title] AND 5.8S[Title] AND 26S[Title] AND \"internal transcribed spacer 2\"[Title]", groups = "plants"),
  metazoa_its = list(term = "Metazoa[Organism] AND 5.8S[Title] AND 28S[Title] AND 18S[Title]", groups = c("vertebrates", "invertebrates")),
  metazoa_18s = list(term = "Metazoa[Organism] AND 18S ribosomal RNA[Title] AND complete sequence[Title] AND 1700:2600[SLEN]", groups = c("vertebrates", "invertebrates")),
  plant_18s = list(term = "Viridiplantae[Organism] AND 18S ribosomal RNA[Title] AND complete sequence[Title] AND 1700:2600[SLEN]", groups = "plants")
)
run_bin <- function(cmd, args) system2(file.path(bin, cmd), args, stdout = TRUE, stderr = FALSE)

tmpl_rds <- file.path(out_dir, "templates.rds")
if (!file.exists(tmpl_rds)) {
  tmpl <- list()
  for (k in names(kinds)) {
    acc <- esearch_acc(kinds[[k]]$term)
    f <- tempfile()
    writeLines(acc, f)
    info <- run_bin("blastdbcmd", c("-db", db, "-entry_batch", f, "-outfmt", shQuote("%a\t%T\t%l")))
    info <- utils::read.delim(text = info, header = FALSE, col.names = c("acc", "taxid", "len"), colClasses = "character")
    info <- info[info$acc %in% acc & !duplicated(info$acc), ]
    tx <- suppressWarnings(get_tax(unique(info$taxid), taxonomy_db = tax_db))
    m <- match(info$taxid, tx$taxid)
    info$species <- tx$species[m]
    info$genus <- tx$genus[m]
    info$family <- tx$family[m]
    lin <- tx$lineage_taxids[m]
    has <- function(id) grepl(paste0("(^|; )", id, "(;|$)"), lin)
    info$group <- ifelse(has(7742), "vertebrates", ifelse(has(33208), "invertebrates", ifelse(has(33090), "plants", NA)))
    info <- info[tx$rank[m] %in% "species" & !is.na(info$genus) & info$group %in% kinds[[k]]$groups &
      !grepl(open_nomenclature_pattern, info$species, perl = TRUE), ]
    info <- info[sample(nrow(info)), ]
    info <- info[!duplicated(info$genus), ]
    info <- utils::head(info, n_templates)
    info$kind <- k
    message(format(Sys.time()), " ", k, ": ", length(acc), " accessions, ", nrow(info), " templates")
    tmpl[[k]] <- info
  }
  tmpl <- do.call(rbind, tmpl)
  # sequences from the local core_nt
  f <- tempfile()
  writeLines(unique(tmpl$acc), f)
  fa <- file.path(out_dir, "templates.fasta")
  run_bin("blastdbcmd", c("-db", db, "-entry_batch", f, "-outfmt", "%f", "-out", fa))
  # every accession of the template entries (identical sequences are merged)
  neg <- run_bin("blastdbcmd", c("-db", db, "-entry_batch", f, "-outfmt", "%a"))
  writeLines(unique(c(tmpl$acc, neg)), file.path(out_dir, "negative.txt"))
  st <- system2(file.path(bin, "blastdb_aliastool"), c(
    "-seqid_file_in", file.path(out_dir, "negative.txt"),
    "-seqid_file_out", file.path(out_dir, "negative.bsl")
  ))
  stopifnot(st == 0L)
  saveRDS(tmpl, tmpl_rds)
}
tmpl <- readRDS(tmpl_rds)
seqs <- readDNAStringSet(file.path(out_dir, "templates.fasta"))
names(seqs) <- sub(" .*", "", names(seqs))
seqs <- seqs[!duplicated(names(seqs))]
message("templates: ", length(seqs), " sequences, ", round(sum(width(seqs)) / 1e6, 1), " Mb")

# 3. in-silico PCR ----------------------------------------------------------
mm <- function(s) if (nchar(s) < 20L) 2L else 3L
pcr <- function(fwd, rev, templ, min_len = 25L, max_len = 1000L) {
  f <- DNAString(fwd)
  rrc <- reverseComplement(DNAString(rev))
  best <- list()
  for (strand in c("+", "-")) {
    s <- if (strand == "+") templ else reverseComplement(templ)
    mf <- vmatchPattern(f, s, max.mismatch = mm(fwd), fixed = c(FALSE, TRUE))
    mr <- vmatchPattern(rrc, s, max.mismatch = mm(rev), fixed = c(FALSE, TRUE))
    for (i in which(lengths(mf) > 0L & lengths(mr) > 0L)) {
      fe <- end(mf[[i]])
      rs <- start(mr[[i]])
      ins <- outer(rs, fe, "-") - 1L # insert between fwd end and rev start
      ok <- which(ins >= min_len & ins <= max_len, arr.ind = TRUE)
      if (!nrow(ok)) next
      j <- ok[which.min(ins[ok]), ]
      len <- ins[j[1L], j[2L]]
      nm <- names(s)[i]
      if (!is.null(best[[nm]]) && nchar(best[[nm]]) <= len) next
      best[[nm]] <- as.character(subseq(s[[i]], start = fe[j[2L]] + 1L, width = len))
    }
  }
  best
}

q_rds <- file.path(out_dir, "queries.rds")
qs <- list()
for (i in seq_len(nrow(p))) {
  for (k in strsplit(p$kinds[i], ";", fixed = TRUE)[[1]]) {
    tk <- tmpl[tmpl$kind == k, ]
    amp <- pcr(p$fwd_seq[i], p$rev_seq[i], seqs[names(seqs) %in% tk$acc], max_len = p$max_len[i])
    if (!length(amp)) next
    d <- data.frame(acc = names(amp), seq = unlist(amp, use.names = FALSE))
    d <- d[!grepl("[^ACGT]", d$seq), ]
    d <- cbind(primer = p$primer[i], marker = p$marker[i], kind = k, d, tk[match(d$acc, tk$acc), c("taxid", "species", "genus", "family", "group")])
    qs[[length(qs) + 1L]] <- d
  }
  message(format(Sys.time()), " ", p$primer[i], ": ", sum(vapply(qs, function(d) sum(d$primer == p$primer[i]), integer(1L))), " amplicons")
}
q <- do.call(rbind, qs)
q <- q[sample(nrow(q)), ]
q <- q[stats::ave(seq_len(nrow(q)), paste(q$primer, q$group), FUN = seq_along) <= n_per_cell, ]
q <- q[order(q$primer, q$group), ]
saveRDS(q, q_rds)
message("queries: ", nrow(q), " (", length(unique(q$seq)), " distinct sequences)")
print(table(q$primer, q$group))
