
blast_cmd <- function(
  query_str,
  db_path,
  num_alignments = "4",
  num_threads = "1",
  blast_type = "blastn",
  perc_id = "80",
  perc_qcov_hsp = "80",
  mt_mode = "2",
  verbose = "silent",
  env_name = "blastr-blast-env"
) {
  query_path <- withr::local_tempfile(
    pattern = "blast_input_",
    fileext = "fasta"
  )
  base::cat(query_str, file = query_path)
  withr::local_envvar(
    .new = list(
      BLAST_USAGE_REPORT = "false",
      NCBI_DONT_USE_NCBIRC = "true",
      NCBI_DONT_USE_LOCAL_CONFIG = "true"
    ),
    action = "replace"
  )
  blast_res <- condathis::run_bin(
    blast_type,
    "-db",
    db_path,
    "-query",
    query_path,
    "-outfmt",
    "6 std qcovhsp staxid stitle qseq",
    "-max_hsps",
    "1",
    "-perc_identity",
    perc_id,
    "-qcov_hsp_perc",
    perc_qcov_hsp,
    "-num_alignments",
    num_alignments,
    "-num_threads",
    num_threads,
    "-mt_mode",
    mt_mode,
    env_name = env_name,
    verbose = verbose,
    error = "continue"
  )
  return(blast_res)
}