library(BLASTr)

# set the number of availble threads to be used (exemplified by the total number of available threads - 2)
# options(BLASTr.num_threads = length(future::availableWorkers()) - 2)

# set the database path (exemplified by the mock blast DB used to be searched with the test ASVs below)
# options(BLASTr.db_path = paste0(fs::path_wd(),"/dev/minimal_db/shortest_minimal_db_BLASTr.fasta"))
# BLASTr.db_path <- paste0(fs::path_wd(),"/dev/minimal_db/shortest_minimal_db_BLASTr.fasta")
# here are 8 ASVs to be tested with the mock blast DB

ASVs_test <- c(
  "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
  "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGTGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
  "GCCAAATTTGTGTTTTGTCCTTCGTTTTTAGTTAATTGTTACTGGCAAATGACTAACGACAAATGATAAATTACTAATAC",
  "AACATTGTATTTTGTCTTTGGGGCCTGGGCAGGTGCAGTAGGAACTTCACTTAGAATAATTATTCGTACTGAGCTTGGGCATCCAGGAAGACTTATCGGGGATGATCAAATCTATAATGTAATTGTTACAGCACATGCATTTGTGATAATTTTTTTTATAGTAATACCTATTATGATT",
  "ACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATG",
  "TTAGCCATAAACATAAAAGTTCACATAACAAGAACTTTTGCCCGAGAACTACTAGCAACAGCTTAAAACTCAAAGGACTTGGCGGTGCTTTATATCCAC"
)



# blaste_res_1ASV <- get_blast_results(asv = "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",num_thread = 10)


future::availableCores()

# options(BLASTr.dbapth = "/data/databases/nt/nt")

blast_res <- BLASTr::parallel_blast(
  asvs = ASVs_test,
  db_path = "/data/databases/core_nt/core_nt",
  # db_path = BLASTr.db_path,
  out_file = NA,
  out_RDS = NA,
  # blast_cmd = "blastn -db {db_path} -outfmt '6 std qcovhsp' -max_hsps 1 -perc_identity 80 -qcov_hsp_perc 80 -num_threads {as.character(num_thread)} -num_alignments {as.character(num_alignments)}",
  total_cores = 77,
  perc_id = 80,
  num_threads = 1,
  perc_qcov_hsp = 80,
  num_alignments = 4,
  blast_type = "blastn"
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
