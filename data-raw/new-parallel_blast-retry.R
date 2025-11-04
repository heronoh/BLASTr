# devtools::load_all()

devtools::load_all()

dna_fasta_path <- fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta")
temp_db_path <- fs::path_temp("minimal_db_blast")
db_path <- temp_db_path
make_blast_db(fasta_path = dna_fasta_path, db_path = temp_db_path)

asvs_string <- "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAG
  AGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC"

run_blast(
  asv = asvs_string,
  db_path = temp_db_path,
  verbose = "full",
  num_threads = 1
)

get_blast_results(asv = asvs_string, db_path = temp_db_path)

tictoc::tic("Parallel BLAST 4 processes with 1 thread each")
parallel_blast(
  asv = asvs_string,
  db_path = temp_db_path,
  total_cores = 4L,
  num_threads = 1L,
  verbose = "full"
)
tictoc::toc()

# Using `blast_cmd`

devtools::load_all()

run_blast(
  asv = asvs_string,
  db_path = temp_db_path,
  verbose = "full",
  num_threads = 1
)

blast_cmd(
  query_str = asvs_string,
  db_path = temp_db_path,
  verbose = "full"
  # strategy = "run_bin"
)

blast_cmd(
  query_str = asvs_string,
  db_path = temp_db_path,
  verbose = "full"
  # strategy = "run"
)


# Benchmarking - microbenchmark
devtools::load_all()
bench_res <- microbenchmark::microbenchmark(
  times = 50,
  check = FALSE,
  run = run_blast(
    asv = asvs_string,
    db_path = temp_db_path,
    verbose = "silent",
    num_threads = 1
  ),
  cmd_run_bin = blast_cmd(
    query_str = asvs_string,
    db_path = temp_db_path,
    verbose = "silent"
    # strategy = "run_bin"
  ),
  cmd_run = blast_cmd(
    query_str = asvs_string,
    db_path = temp_db_path,
    verbose = "silent"
    # strategy = "run"
  )
)

bench_res

plot(bench_res)

# Benchmarking - bench::mark
# Single string input 1 core
devtools::load_all()
mark_res <- bench::mark(
  iterations = 100,
  check = FALSE,
  run = run_blast(
    asv = asvs_string,
    db_path = temp_db_path,
    verbose = "silent",
    num_threads = 1
  ),
  cmd_run_bin = blast_cmd(
    query_str = asvs_string,
    db_path = temp_db_path,
    verbose = "silent"
    # strategy = "run_bin"
  ),
  cmd_run = blast_cmd(
    query_str = asvs_string,
    db_path = temp_db_path,
    verbose = "silent"
    # strategy = "run"
  )
)

print(mark_res)
plot(mark_res)

# Test CTRL-C interruption handling in parallel_blast2

devtools::load_all()

mirai::status()
tictoc::tic("Parallel BLAST 1000 processes with 8 processes and 1 thread each - new")
par_res_1000_new <- parallel_blast2(
  query_seqs = rep(asvs_string, 1000L),
  db_path = temp_db_path,
  total_cores = 8L,
  num_threads = 1L,
  verbose = "full"
)
tictoc::toc()

tictoc::tic("Parallel BLAST 1000 processes with 8 processes and 1 thread each - old")
par_res_1000_old <- parallel_blast(
  asvs = rep(asvs_string, 1000L),
  db_path = temp_db_path,
  total_cores = 8L,
  num_threads = 1L,
  verbose = "full"
)
tictoc::toc()

# =============================================================================
# New straight implementation of parallel_blast
# =============================================================================
devtools::load_all()
dna_fasta_path <- fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta")
temp_db_path <- fs::path_temp("minimal_db_blast")
make_blast_db(fasta_path = dna_fasta_path, db_path = temp_db_path)

asvs_string <- "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAG
  AGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC"

db_path <- temp_db_path
db_files_info <- fs::dir_info(
  path = fs::path_dir(db_path),
  regexp = fs::path_file(db_path)
)[c("path", "size")]
# str(db_files_info)
db_hash <- rlang::hash(db_files_info)
input_hash <- rlang::hash(asvs_string)
query_seq <- stringr::str_replace_all(asvs_string, "\\s", "")
query_hash <- rlang::hash(query_seq)
tibble::tibble(
  # db_file = db_files,
  db_hash = db_hash,
  input_hash = input_hash,
  query_hash = query_hash
)

tictoc::tic("Parallel BLAST 4 processes with 1 thread each - new")
par_res_10_new <- parallel_blast(
  query_seqs = rep(asvs_string, 10L),
  db_path = temp_db_path,
  total_cores = 4L,
  num_threads = 1L,
  verbose = "silent"
)
tictoc::toc()

tictoc::tic("Parallel BLAST 4 processes with 1 thread each - old")
par_res_10_old <- parallel_blast_old(
  asvs = rep(asvs_string, 10L),
  db_path = temp_db_path,
  total_cores = 4L,
  num_threads = 1L,
  verbose = "silent"
)
tictoc::toc()

waldo::compare(par_res_10_old, par_res_10_new, max_diffs = 100)
colnames(par_res_10_old)
colnames(par_res_10_new)
par_res_10_old$`1_subject header`
par_res_10_new$`1_subject header`

mark_par_res <- bench::mark(
  iterations = 10,
  check = TRUE,
  old = parallel_blast_old(
    asvs = rep(asvs_string, 10L),
    db_path = temp_db_path,
    total_cores = 4L,
    num_threads = 1L,
    verbose = "silent"
  ),
  new = parallel_blast(
    query_seqs = rep(asvs_string, 10L),
    db_path = temp_db_path,
    total_cores = 4L,
    num_threads = 1L,
    verbose = "silent"
  )
)

print(mark_par_res)
plot(mark_par_res)


fail_query_seq <- "AGAGTACTACAAGTGCTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC"


blast_cmd(
  query_str = fail_query_seq,
  db_path = temp_db_path,
  num_alignments = 10L,
  verbose = "silent"
)

parallel_blast(
  query_seqs = fail_query_seq,
  db_path = temp_db_path,
  total_cores = 2L,
  num_threads = 1L,
  num_alignments = 10L,
  verbose = "full"
)

res_df <- parallel_blast(
  query_seqs = asvs_string,
  db_path = temp_db_path,
  total_cores = 2L,
  num_threads = 1L,
  num_alignments = 10L,
  verbose = "full"
)

# View(res_df)

par_res <- parallel_blast(
  query_seqs = rep(asvs_string, 1000L),
  db_path = temp_db_path,
  total_cores = 4L,
  num_threads = 1L,
  num_alignments = 10L
)

# Generate random combinations of letters A, T, C, G

random_seqs <- function(num = 100L) {
  set.seed(123)
  base::sapply(base::seq_len(num), function(x) {
    paste0(sample(c("A", "T", "C", "G"), size = 50, replace = TRUE), collapse = "")
  })
}


no_match_query_seq <- "AGAGTACTACAAGTGCTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC"
error_header_seq <- ">XAVADADAA"
error_char_seq <- "ZZZaaaa1345!!!@@@"

par_res <- parallel_blast(
  query_seqs = rep(c(asvs_string, no_match_query_seq, error_header_seq, error_char_seq), 2L),
  db_path = temp_db_path,
  total_cores = 1L,
  num_threads = 1L,
  num_alignments = 10L,
  retry_times = 3L,
  verbose = "progress"
)

par_res <- parallel_blast(
  query_seqs = random_seqs(10000L),
  db_path = temp_db_path,
  total_cores = 1L,
  num_threads = 1L,
  num_alignments = 10L,
  retry_times = 3L,
  verbose = "progress"
)

par_res <- parallel_blast(
  query_seqs = random_seqs(100L),
  db_path = temp_db_path,
  total_cores = 1L,
  num_threads = 1L,
  num_alignments = 10L,
  retry_times = 3L,
  verbose = "output"
)

par_res <- parallel_blast(
  query_seqs = random_seqs(100L),
  db_path = temp_db_path,
  total_cores = 1L,
  num_threads = 1L,
  num_alignments = 10L,
  retry_times = 3L,
  verbose = "output"
)

par_res <- parallel_blast(
  query_seqs = random_seqs(100L),
  db_path = temp_db_path,
  total_cores = 1L,
  num_threads = 1L,
  num_alignments = 10L,
  retry_times = 3L,
  verbose = "cmd"
)


attributes(par_res)

# This should error
par_old_res <- parallel_blast_old(
  asvs = rep(c(asvs_string, no_match_query_seq, error_header_seq, error_char_seq), 2L),
  db_path = temp_db_path,
  total_cores = 4L,
  num_threads = 1L,
  num_alignments = 10L,
  verbose = "full"
)

blast_cmd(
  query_str = no_match_query_seq,
  db_path = temp_db_path
)
blast_cmd(
  query_str = error_header_seq,
  db_path = temp_db_path
)
blast_cmd(
  query_str = error_char_seq,
  db_path = temp_db_path
)
