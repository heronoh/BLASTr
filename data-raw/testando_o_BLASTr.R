# devtools::install_github("heronoh/BLASTr")
library(BLASTr)

devtools::load_all()

install_dependencies(force = TRUE, verbose = "full")

packageVersion("BLASTr")

packageVersion("condathis")

packageVersion("purrr")

# set the number of availble threads to be used (exemplified by the total number of available threads - 2)
# options(BLASTr.num_threads = length(future::availableWorkers()) - 2)

# set the database path (exemplified by the mock blast DB used to be searched with the test ASVs below)
# options(BLASTr.db_path = paste0(fs::path_wd(), "/data-raw/minimal_db/shortest_minimal_db_BLASTr.fasta"))

BLASTr.db_path <- paste0(fs::path_wd(), "/data-raw/minimal_db/shortest_minimal_db_BLASTr.fasta")
# here are 8 ASVs to be tested with the mock blast DB
# options(BLASTr.dbapth = "/data/databases/nt/nt")

# fs::path_package()



ASVs_test <- c(
  "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
  "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGTGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
  "GCCAAATTTGTGTTTTGTCCTTCGTTTTTAGTTAATTGTTACTGGCAAATGACTAACGACAAATGATAAATTACTAATAC",
  "AACATTGTATTTTGTCTTTGGGGCCTGGGCAGGTGCAGTAGGAACTTCACTTAGAATAATTATTCGTACTGAGCTTGGGCATCCAGGAAGACTTATCGGGGATGATCAAATCTATAATGTAATTGTTACAGCACATGCATTTGTGATAATTTTTTTTATAGTAATACCTATTATGATT",
  "XXXXAAANNN",
  "ACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATG",
  "TTAGCCATAAACATAAAAGTTCACATAACAAGAACTTTTGCCCGAGAACTACTAGCAACAGCTTAAAACTCAAAGGACTTGGCGGTGCTTTATATCCAC"
)

# blaste_res_ASV <- get_blast_results(asv = "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",num_thread = 10)

parallelly::availableCores()
# New mirai based parallel
mirai::daemons(n = 4L, .compute = "blastr-cpu")
mirai::require_daemons(.compute = "blastr-cpu")

mirai::status(.compute = "blastr-cpu")
mirai::status(.compute = "blastr-cpu")$connections

purrr::in_parallel()
purrr::map()
purrr:::map_()

# So for `parallel_blast` case

# Test single process
blast_single_res <- parallel_blast(
  asvs = ASVs_test,
  # db_path = "/data/databases/nt/nt",
  db_path = BLASTr.db_path,
  out_file = NA,
  out_RDS = NA,
  total_cores = 1L,
  perc_id = 80L,
  num_threads = 1L,
  perc_qcov_hsp = 80L,
  num_alignments = 4L,
  blast_type = "blastn",
  verbose = "full"
)

# Testing multiple processes
blast_multi_res <- parallel_blast(
  asvs = ASVs_test,
  # db_path = "/data/databases/nt/nt",
  db_path = BLASTr.db_path,
  out_file = NA,
  out_RDS = NA,
  total_cores = 2L,
  perc_id = 80,
  num_threads = 1,
  perc_qcov_hsp = 80,
  num_alignments = 4,
  blast_type = "blastn",
  verbose = "full"
)

identical(blast_multi_res, blast_single_res)

teste <- readRDS(file = "/home/heron/prjcts/omics/metaseqs/blast_out.RDS")

# Close connections
mirai::daemons(n = 0L, .compute = "blastr-cpu")


# benchmarking

## To Run Background Job

# rstudioapi::jobRunScript(path = "heron-otu-blastn.R", workingDir = ".", importEnv = TRUE, exportEnv = "R_GlobalEnv")

# ----
