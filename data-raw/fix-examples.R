# As of v0.1.5 - 2025-10-28
# + Most examples do not really run because the database path is not available.
# Fixed on v0.1.6

devtools::load_all()

dna_fasta_path <- fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta")
temp_db_path <- fs::path_temp("minimal_db_blast")
make_blast_db(fasta_path = dna_fasta_path, db_path = temp_db_path)

asvs_string <- "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAG
  AGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC"

# This example code fails if db is not build yet
blast_res <- run_blast(
  asv = asvs_string,
  db_path = dna_fasta_path,
  verbose = "full"
)


run_blast(
  asv = "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAG
  AGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
  db_path = temp_db_path,
  verbose = "silent",
  num_thread = 0
)
# Using temp formatted blast db
blast_res <- run_blast(
  asv = "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAG
  AGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
  db_path = temp_db_path,
  perc_id = 80,
  num_threads = 1,
  perc_qcov_hsp = 80,
  num_alignments = 2,
  verbose = FALSE
)

blast_res <- run_blast(asv = asvs_string, db_path = temp_db_path)
blast_res


get_blast_results(asv = asvs_string, db_path = temp_db_path)

get_fasta_header()


parallel_blast(
  asv = asvs_string,
  db_path = temp_db_path,
  verbose = "full"
)

get_fasta_header()

get_fasta_header(id = "AP011979.1", db_path = temp_db_path)
