library(BLASTr)


# condathis::create_env("bioconda::entrez-direct==22.4", env_name = "entrez-env")

# condathis::create_env("bioconda::blast==2.16", env_name = "blast-env")

install_dependencies()

packageVersion("BLASTr")
# if (isFALSE(nzchar(cmd_bin))) {
#
#
#   cmd_bin <- Sys.which(cmd)
#   if (isFALSE(nzchar(cmd_bin))) {
#
#     if (!condathis::env_exists("blast-env")) {
#       condathis::create_env("bioconda::blast==2.16", env_name = "blast-env")
#     }
# }

# set the number of availble threads to be used (exemplified by the total number of available threads - 2)
# options(BLASTr.num_threads = length(future::availableWorkers()) - 2)

# set the database path (exemplified by the mock blast DB used to be searched with the test ASVs below)
# options(BLASTr.db_path = paste0(fs::path_wd(), "/data-raw/minimal_db/shortest_minimal_db_BLASTr.fasta"))

BLASTr.db_path <- paste0(fs::path_wd(), "/data-raw/minimal_db/shortest_minimal_db_BLASTr.fasta")
# here are 8 ASVs to be tested with the mock blast DB

ASVs_test <- c(
  "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
  "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGTGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
  "GCCAAATTTGTGTTTTGTCCTTCGTTTTTAGTTAATTGTTACTGGCAAATGACTAACGACAAATGATAAATTACTAATAC",
  "AACATTGTATTTTGTCTTTGGGGCCTGGGCAGGTGCAGTAGGAACTTCACTTAGAATAATTATTCGTACTGAGCTTGGGCATCCAGGAAGACTTATCGGGGATGATCAAATCTATAATGTAATTGTTACAGCACATGCATTTGTGATAATTTTTTTTATAGTAATACCTATTATGATT",
  "ACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATG",
  "TTAGCCATAAACATAAAAGTTCACATAACAAGAACTTTTGCCCGAGAACTACTAGCAACAGCTTAAAACTCAAAGGACTTGGCGGTGCTTTATATCCAC"
)



# blaste_res_ASV <- get_blast_results(asv = "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",num_thread = 10)


future::availableCores()

# options(BLASTr.dbapth = "/data/databases/nt/nt")

blast_res <- parallel_blast(
  asvs = ASVs_test,
  # db_path = "/data/databases/nt/nt",
  db_path = BLASTr.db_path,
  out_file = NA,
  out_RDS = NA,
  total_cores = 4,
  perc_id = 80,
  num_threads = 1,
  perc_qcov_hsp = 80,
  num_alignments = 4,
  blast_type = "blastn",
  verbose = "full"
)

# blast_res1 <- BLASTr::run_blast(
#   asv = ASVs_test[3],
#   db_path = "/data/databases/nt/nt",
#   perc_id = 80,
#   num_thread = 1,
#   perc_qcov_hsp = 80,
#   num_alignments = 2
#   # ,
#   # blast_type = "blastn"
# )

teste <- readRDS(file = "/home/heron/prjcts/omics/metaseqs/blast_out.RDS")

# benchmarking

## To Run Background Job

# rstudioapi::jobRunScript(path = "heron-otu-blastn.R", workingDir = ".", importEnv = TRUE, exportEnv = "R_GlobalEnv")

# ----
