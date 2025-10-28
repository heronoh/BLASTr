# As of v0.1.5 - 2025-10-28
# + Most examples do not really run because the database path is not available.

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

# Single Blast process with 10 Threads
tictoc::tic("Single BLAST process with 10 threads")
parallel_blast(
  asv = asvs_string,
  db_path = temp_db_path,
  total_cores = 1L,
  num_threads = 10L,
  verbose = "full"
)
tictoc::toc()


parallel_blast(
  asv = asvs_string,
  db_path = temp_db_path,
  total_cores = 1L,
  num_thread = 10L,
  verbose = "silent"
)

# 10 Blast processes with 2 Threads each
parallel_blast(
  asv = asvs_string,
  db_path = temp_db_path,
  total_cores = 10L,
  num_threads = 2L,
  verbose = "full"
)

# 10 Blast processes with 2 Threads each
tictoc::tic("Parallel BLAST 10 processes with 2 threads each")
parallel_blast(
  asv = rep(asvs_string, 10),
  db_path = temp_db_path,
  total_cores = 10L,
  num_threads = 2L,
  verbose = "output"
)
tictoc::toc()
