# As of v0.1.5 - 2025-10-28
# + `blastn` is returning warning with exit code 3 and messages are not properly captured.


devtools::load_all()

dna_fasta_path <- fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta")
temp_db_path <- fs::path_temp("minimal_db_blast")
make_blast_db(fasta_path = dna_fasta_path, db_path = temp_db_path)

asvs_string <- "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAG
  AGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC"

# This example code fails if db is not build yet
# + but message is not that informative.
# + Error and Warning are globbed in a single line.
blast_res <- run_blast(
  asv = asvs_string,
  db_path = dna_fasta_path,
  verbose = "full"
)

# Warning is still captured at stderr redirection
run_blast(
  asv = "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAG
  AGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
  db_path = temp_db_path,
  verbose = "silent"
)
