# Objetivo: dar taxid ao banco BLAST minimo do pacote
# (data-raw/minimal_db/shortest_minimal_db_BLASTr.fasta), para que
# parallel_blast()/run_blast() retornem `staxid` tambem no banco de teste.
#
# Fluxo:
#   1. ler o FASTA minimo
#   2. buscar cada sequencia no core_nt para descobrir o taxid do melhor hit
#   3. escrever o taxid map (lcl|<ID> <taxid>)
#   4. recriar o banco minimo com -taxid_map
#   5. validar: BLAST contra o banco minimo deve trazer staxid
#
# Comandos de referencia (shell):
# echo -e 'GAGGCCCCGCTACGAAAGTAACTTTATCACCTCCGACCCCACGAAAGCTAAGAAACAAACTGGGATTAGATACCCCACTATGCTTAGCCCTAAACCTAGATATTTATTCTACAAAAAATATCCGCCAGGGAACTACGAGCGCTAGCTTAAAACCCAAAGGACTTGACGGTGTCTCAGACCCACCTAGAGGAGCCTGTTCTAGAACCGATAATCCACGTTTAACCTCACCACCCCTTGTTTTCCCCGCCTATATACCGCCGTCGCAAGCTTACCCTGTGA' | blastn -db "$BLASTR_LOCAL_DB" -outfmt "6 std qcovhsp ssciname" -max_hsps 3 -perc_identity 80 -num_threads 80 -num_alignments 3
# echo -e 'GAGGCCCCGCTACGAAAGTAACTTTATCACCTCCGACCCCACGAAAGCTAAGAAACAAACTGGGATTAGATACCCCACTATGCTTAGCCCTAAACCTAGATATTTATTCTACAAAAAATATCCGCCAGGGAACTACGAGCGCTAGCTTAAAACCCAAAGGACTTGACGGTGTCTCAGACCCACCTAGAGGAGCCTGTTCTAGAACCGATAATCCACGTTTAACCTCACCACCCCTTGTTTTCCCCGCCTATATACCGCCGTCGCAAGCTTACCCTGTGA' | blastn -db "$BLASTR_CORE_NT_DB" -outfmt "6 std qcovhsp ssciname" -max_hsps 3 -perc_identity 80 -num_threads 80 -num_alignments 3

devtools::load_all(path = ".")

core_nt_path <- Sys.getenv("BLASTR_CORE_NT_DB")
minimal_db_fasta <- "data-raw/minimal_db/shortest_minimal_db_BLASTr.fasta"
taxid_map_path <- "data-raw/minimal_db/shortest_minimal_db_BLASTr.txt"

# O BLAST precisa achar o taxdb (taxdb.btd/.bti) do core_nt.
# Equivale a: export BLASTDB="$(dirname "$BLASTR_CORE_NT_DB"):$BLASTDB"
Sys.setenv(
  BLASTDB = paste(
    c(dirname(core_nt_path), Sys.getenv("BLASTDB")),
    collapse = ":"
  ) |>
    sub(pattern = ":$", replacement = "")
)

# 1. Ler o FASTA minimo ------------------------------------------------------

ref_seqs <- Biostrings::readDNAStringSet(filepath = minimal_db_fasta)

seqs_tbl <- tibble::tibble(
  "names" = stringr::str_replace_all(names(ref_seqs), pattern = "\\.", replacement = "_"),
  "Sequence" = as.character(ref_seqs)
)

# nomes nao podem repetir (viram IDs no makeblastdb -parse_seqids)
stopifnot(!any(duplicated(seqs_tbl$names)))

# sequencias usadas como query nos passos seguintes
seqs <- as.character(ref_seqs) |> unname()

# 2. Descobrir o taxid de cada sequencia no core_nt --------------------------

blast_res_test <- parallel_blast(
  db_path = core_nt_path,
  query_seqs = seqs,
  total_cores = 60,
  perc_id = 80,
  num_threads = 2,
  perc_qcov_hsp = 80,
  num_alignments = 3,
  blast_type = "blastn"
)

# 3. Escrever o taxid map ----------------------------------------------------

seqs_tbl <- seqs_tbl |>
  dplyr::left_join(
    blast_res_test,
    by = "Sequence",
    relationship = "one-to-one"
  ) |>
  dplyr::mutate("IDs" = stringr::str_remove_all(names, pattern = " .*$"))

paste0("lcl|", seqs_tbl$IDs, " ", seqs_tbl$`1_staxid`) |>
  write(taxid_map_path)

# 4. Recriar o banco minimo com taxid ----------------------------------------
# Equivale a:
# makeblastdb -in shortest_minimal_db_BLASTr.fasta -parse_seqids -blastdb_version 5 -taxid_map shortest_minimal_db_BLASTr.txt -dbtype nucl -hash_index

make_blast_db(
  fasta_path = minimal_db_fasta,
  verbose = "full",
  db_type = "nucl",
  taxid_map = taxid_map_path
)

# 5. Validar -----------------------------------------------------------------
# Referencia em shell:
# echo -e 'AACATTGTATTTTGTCTTTGGGGCCTGGGCAGGTGCAGTAGGAACTTCACTTAGAATAATTATTCGTACTGAGCTTGGGCATCCAGGAAGACTTATCGGGGATGATCAAATCTATAATGTAATTGTTACAGCACATGCATTTGTGATAATTTTTTTTATAGTAATACCTATTATGATT' | blastn -db 'data-raw/minimal_db/shortest_minimal_db_BLASTr.fasta' -outfmt "6 std qcovhsp ssciname" -max_hsps 3 -perc_identity 80 -num_threads 80 -num_alignments 3

blast_results <- parallel_blast(
  query_seqs = seqs,
  db_path = minimal_db_fasta,
  # db_path = core_nt_path, # para comparar com o core_nt
  num_alignments = 3,
  total_cores = 2 # Number of BLAST processes to use
)

blast_results

# todas as sequencias deveriam ter staxid no banco minimo
stopifnot(!anyNA(blast_results$`1_staxid`))

# run_blast() escreve `query_seqs` com cat(): um vetor de varias sequencias
# vira UMA query concatenada por espacos. Para buscar mais de uma, passar
# um FASTA ja formatado em uma unica string:
query_fasta <- paste0(">q", 3:4, "\n", seqs[3:4], collapse = "\n")

run_blast(
  query_seqs = query_fasta,
  db_path = minimal_db_fasta,
  # db_path = core_nt_path,
  num_alignments = 3
)
