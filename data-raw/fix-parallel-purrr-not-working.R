# As of v0.1.5 - 2025-10-28
# + parallel blast is not really parallelizing anything
# + probably because of how mirai::daemons() are set inside the function.

# Fixed on v0.1.6

# devtools::load_all()
# library(BLASTr)

run_example_parallel <- function(total_cores, num_threads = 1L) {
  # library(BLASTr)
  devtools::load_all()
  dna_fasta_path <- fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta")
  temp_db_path <- withr::local_tempfile()
  make_blast_db(fasta_path = dna_fasta_path, db_path = temp_db_path)
  # query_string <- "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCA"
  query_string <- "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAG
AGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC"


  res <- parallel_blast(
    asv = rep(query_string, 1000),
    db_path = temp_db_path,
    total_cores = total_cores,
    num_threads = num_threads,
    verbose = "silent"
  )
  # withr::deferred_clear()
  return(res)
}
# TODO: @luciorq Test again after implementing the argument mt_mode
tictoc::tic("Example Parallel BLAST with 20 Cores and 1 Thread each")
run_example_parallel(20L, 1L)
tictoc::toc()

tictoc::tic("Example Parallel BLAST with 4 Cores and 1 Thread each")
run_example_parallel(4L, 1L)
tictoc::toc()

tictoc::tic("Example Parallel BLAST with 4 Cores and 1 Thread each")
run_example_parallel(4L, 4L)
tictoc::toc()



dna_fasta_path <- fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta")
temp_db_path <- fs::path_temp("minimal_db_blast")
fs::dir_exists(temp_db_path)
make_blast_db(fasta_path = dna_fasta_path, db_path = temp_db_path)

asvs_string <- "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAG
AGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC"


# This should Error with a fail for partial matching
parallel_blast(
  asv = rep(asvs_string, 100),
  db_path = temp_db_path,
  total_cores = 20L,
  num_threads = 20L,
  verbose = "silent"
)


# Tests with injected sleep of 3 seconds

# Single Blast process with 10 Threads - Should be around 3 seconds
# + 3.209 secs
tictoc::tic("Single BLAST process with 10 threads")
parallel_blast(
  asv = asvs_string,
  db_path = temp_db_path,
  total_cores = 1L,
  num_threads = 10L,
  verbose = "full"
)
tictoc::toc()


# 10 Blast processes with 2 Threads each - around 3 secs as it is just 1 sequence
# + 3.48 secs - Overhead of creating daemons - Should be fixed in v0.1.6
tictoc::tic("Parallel BLAST 10 processes with 2 threads each")
parallel_blast(
  asv = asvs_string,
  db_path = temp_db_path,
  total_cores = 10L,
  num_threads = 2L,
  verbose = "full"
)
tictoc::toc()

# 10 Blast processes with 2 Threads each
# +  31.516 secs - it is not really parallelizing
tictoc::tic("Parallel BLAST 10 processes with 2 threads each")
parallel_blast(
  asv = rep(asvs_string, 10),
  db_path = temp_db_path,
  total_cores = 10L,
  num_threads = 2L,
  verbose = "output"
)
tictoc::toc()


# 10 Blast processes with 2 Threads each - Manual Mirai configuration
# + 3.881 secs - Should be around 3 seconds

# parallel_strategy <- function(n = 10L) {
#   x <- mirai::daemons(n = n, .compute = NULL)
#   invisible(x)
# }
# # mirai::daemons(n = 10L, .compute = "blastr-cpu")
# # mirai::daemons(n = 10L)

# parallel_strategy(n = 8L)
tictoc::tic("Parallel BLAST 10 processes with 8 threads each - manual Mirai configuration")
res <- parallel_blast(
  asv = rep(asvs_string, 8),
  db_path = temp_db_path,
  total_cores = 10L,
  num_threads = 8L,
  verbose = "output"
)
tictoc::toc()

# mirai::daemons(n = 0L, .compute = "blastr-cpu")
# mirai::status(.compute = "blastr-cpu")

mirai::daemons(0L)
mirai::status()

purrr::in_parallel()

purrr:::map_

purrr:::mmap_
