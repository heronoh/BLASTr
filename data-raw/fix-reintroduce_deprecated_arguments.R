# For v0.1.7
# Reintroduce deprecated arguments to parallel_blast

# Deprecating arguments example
# From <https://design.tidyverse.org/argument-clutter.html#how-do-i-remediate-past-mistakes>
my_fun <- function(x, y, opts = my_fun_opts(), opt1 = lifecycle::deprecated(), opt2 = lifecycle::deprecated()) {
  # browser()
  if (lifecycle::is_present(opt1)) {
    lifecycle::deprecate_warn("1.0.0", "my_fun(opt1)", "my_fun_opts(opt1)")
    opts$opt1 <- opt1
  }
  if (lifecycle::is_present(opt2)) {
    lifecycle::deprecate_warn("1.0.0", "my_fun(opt2)", "my_fun_opts(opt2)")
    opts$opt2 <- opt2
  }
}

# -----
# Bring back deprecated arguments to parallel_blast
devtools::load_all()
dna_fasta_path <- fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta")
temp_db_path <- fs::path_temp("minimal_db_blast")
make_blast_db(fasta_path = dna_fasta_path, db_path = temp_db_path)

asvs_string <- "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAG
  AGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC"

parallel_blast(
  query_seqs = asvs_string,
  db_path = temp_db_path,
  total_cores = 1L,
  num_threads = 1L
)


parallel_blast(
  query_seqs = asvs_string,
  db_path = temp_db_path,
  total_cores = 1L,
  num_threads = 1L,
  out_file = fs::path_temp("blast_results.csv")
)

parallel_blast(
  query_seqs = asvs_string,
  db_path = temp_db_path,
  total_cores = 1L,
  num_threads = 1L,
  # out_file = fs::path_temp("blast_results.csv")
  out_RDS = fs::path_temp("blast_results.rds")
)

# Should fail
parallel_blast(
  asv = asvs_string,
  db_path = temp_db_path,
  total_cores = 1L,
  num_threads = 1L
)

parallel_blast(
  asvs = asvs_string,
  db_path = temp_db_path,
  total_cores = 1L,
  num_threads = 1L
)
