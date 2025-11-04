library(BLASTr)
dna_fasta_path <- fs::path_package(
  "BLASTr", "extdata", "minimal_db_blast",
  ext = "fasta"
)
temp_db_path <- fs::path_temp("minimal_db_blast")
make_blast_db(fasta_path = dna_fasta_path, db_path = temp_db_path)
asvs_test <- readLines(fs::path_package("BLASTr", "extdata", "asvs_test.txt"))

no_match_query_seq <- "AGAGTACTACAAGTGCTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC"
error_header_seq <- ">XAVADADAA"
error_char_seq <- "ZZZaaaa1345!!!@@@"

asvs_to_run <- c(
  asvs_test,
  no_match_query_seq,
  error_header_seq,
  error_char_seq
)

blast_res <- parallel_blast(
  query_seqs = asvs_to_run,
  db_path = temp_db_path,
  total_cores = 1L,
  num_threads = 1L,
  num_alignments = 10L,
  retry_times = 3L
)

print(blast_res)

exit_codes(blast_res)

attributes(blast_res)
