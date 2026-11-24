# Writes inst/extdata/blast_output_example.tsv: the raw tabular output of
# BLAST+ (the -outfmt of parallel_blast()) for the sequences of
# inst/extdata/seqs_test.txt against the database made from
# inst/extdata/minimal_db_blast.fasta, with the search options of
# parallel_blast() (identity 80, coverage 50, 100 subjects, one HSP each).
# The tests read it to check the steps after the search (reading, ranking,
# summary) where BLAST+ is not installed, and compare it with a new search
# where it is.
# Written with BLAST+ 2.17.0. Run from the package root, with BLAST+
# installed (install_dependencies()).
devtools::load_all(".")
fasta <- "inst/extdata/minimal_db_blast.fasta"
db <- fs::path_temp("minimal_db_blast")
make_blast_db(fasta, db_path = db)
seqs <- unique(readLines("inst/extdata/seqs_test.txt"))
query <- fs::path_temp("query", ext = "fasta")
writeLines(paste0(">", sprintf("s%06d", seq_along(seqs)), "\n", seqs), query)
res <- run_blast_program("blastn", c(
  "-db", db, "-query", query, "-outfmt", blast_outfmt(),
  "-max_hsps", "1", "-perc_identity", "80", "-qcov_hsp_perc", "50", "-num_alignments", "100"
))
cat(res$stdout, file = "inst/extdata/blast_output_example.tsv")
