

# echo -e 'GAGGCCCCGCTACGAAAGTAACTTTATCACCTCCGACCCCACGAAAGCTAAGAAACAAACTGGGATTAGATACCCCACTATGCTTAGCCCTAAACCTAGATATTTATTCTACAAAAAATATCCGCCAGGGAACTACGAGCGCTAGCTTAAAACCCAAAGGACTTGACGGTGTCTCAGACCCACCTAGAGGAGCCTGTTCTAGAACCGATAATCCACGTTTAACCTCACCACCCCTTGTTTTCCCCGCCTATATACCGCCGTCGCAAGCTTACCCTGTGA' | blastn -db '/home/gabriel/projetos/db-LGC/old/database/blastr/LGC12Sdb_complete_noGaps-2024-09-20.fasta' -outfmt "6 std qcovhsp ssciname" -max_hsps 3 -perc_identity 80 -num_threads 80 -num_alignments 3
# echo -e 'GAGGCCCCGCTACGAAAGTAACTTTATCACCTCCGACCCCACGAAAGCTAAGAAACAAACTGGGATTAGATACCCCACTATGCTTAGCCCTAAACCTAGATATTTATTCTACAAAAAATATCCGCCAGGGAACTACGAGCGCTAGCTTAAAACCCAAAGGACTTGACGGTGTCTCAGACCCACCTAGAGGAGCCTGTTCTAGAACCGATAATCCACGTTTAACCTCACCACCCCTTGTTTTCCCCGCCTATATACCGCCGTCGCAAGCTTACCCTGTGA' | blastn -db '/data/databases/core_nt/core_nt' -outfmt "6 std qcovhsp ssciname" -max_hsps 3 -perc_identity 80 -num_threads 80 -num_alignments 3
#
#
#


seqs <- Biostrings::readDNAStringSet(filepath = "data-raw/minimal_db/shortest_minimal_db_BLASTr.fasta")


seqs_tbl <- tibble::tibble("names" = stringr::str_replace_all(names(seqs),pattern = "\\.",
                                                              replacement = "_"),
                           "Sequence" = as.character(seqs))

# tem que exportar essa variavel pra funcionar!!!!!
# export BLASTDB=/data/databases/core_nt:/data/databases/core_nt:$BLASTDB


"/data/databases/core_nt/core_nt"

devtools::load_all(path = "/home/noreh/prjcts/BLASTr/R")

#----



# blast_res_test <- parallel_blast(
#   db_path = "/data/databases/core_nt/core_nt",
#   query_seqs = as.character(seqs),
#   # asvs = asvs_blast_all,
#   # asvs = asvs_blast_all_new,
#   # out_file = paste0(blast_path,"/blast_out_res_1.csv"),
#   # out_RDS = paste0(blast_path,"/blast_out_res_1.RDS"),
#   total_cores = 60,
#   perc_id = 80,
#   num_threads = 2,
#   perc_qcov_hsp = 80,
#   num_alignments = 3,
#   blast_type = "blastn"
# )


seqs_tbl$names

seqs_tbl$names |> duplicated()


seqs_tbl <- seqs_tbl |>
  dplyr::left_join(blast_res_test,
                   by = "Sequence",
                   relationship = "one-to-one"
                   ) |>
  dplyr::mutate("IDs" = stringr::str_remove_all(names, pattern = " .*$"))
  # dplyr::mutate("IDs" = names)

paste0("lcl|",seqs_tbl$IDs," ",seqs_tbl$`1_staxid`) |>
  cbind() |>

    write("data-raw/minimal_db/shortest_minimal_db_BLASTr.txt")


# makeblastdb -in shortest_minimal_db_BLASTr.fasta -parse_seqids -blastdb_version 5 -taxid_map shortest_minimal_db_BLASTr.fasta.txt -dbtype nucl -hash_index



# echo -e 'AACATTGTATTTTGTCTTTGGGGCCTGGGCAGGTGCAGTAGGAACTTCACTTAGAATAATTATTCGTACTGAGCTTGGGCATCCAGGAAGACTTATCGGGGATGATCAAATCTATAATGTAATTGTTACAGCACATGCATTTGTGATAATTTTTTTTATAGTAATACCTATTATGATT' | blastn -db '/home/noreh/prjcts/BLASTr/data-raw/minimal_db/shortest_minimal_db_BLASTr.fasta' -outfmt "6 std qcovhsp ssciname" -max_hsps 3 -perc_identity 80 -num_threads 80 -num_alignments 3
#

condathis::run()
# ----

make_blast_db(fasta_path = "data-raw/minimal_db/shortest_minimal_db_BLASTr.fasta",

              verbose = "full",
              db_type = "nucl",
              taxid_map = "data-raw/minimal_db/shortest_minimal_db_BLASTr.txt")



blast_results <- parallel_blast(
  query_seqs = asvs,
  db_path = "data-raw/minimal_db/shortest_minimal_db_BLASTr.fasta",
  num_alignments = 3,
  total_cores = 2 # Number of BLAST processes to use
)

blast_results

