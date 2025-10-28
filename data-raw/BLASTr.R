# # BLAST parser for parallel ASV identification using BLAST suite and local Databases
#
#
#
#
# #1 - Load required functions ----
# #1a - shell_exec ----
#
# ###                          Execute shell commands                          ###
# ## from Lucio RQ - https://rdrr.io/github/luciorq/luciolib/man/shell_exec.html##
# shell_exec <- function(cmd, .envir = parent.frame()) {
#   if (!requireNamespace("processx", quietly = TRUE)) {
#     rlang::abort(message = "Package `processx` package is not installed.")
#   }
#   cmd_res <- processx::run(
#     command = "bash",
#     args = c("-c", glue::glue(cmd, .envir = .envir)), echo_cmd = FALSE
#   )
#   return(cmd_res)
# }
#
# #1b - run_blast ----
# ###                  function to run blast for each ASV/OTU                  ###
# run_blast <- function(asv, db_path = "/data/databases/nt/nt", num_alignments = 3, num_threads = 40) {
#
#   #BLAST command
#   blast_cmd <- "echo -e '>seq1\n{asv}' | blastn -db {db_path} -outfmt '6 std qcovhsp' -max_hsps 1 -perc_identity 80 -qcov_hsp_perc 80 -num_threads {as.character(num_thread)} -num_alignments {as.character(num_alignments)}"
#
#
#   blast_res <- shell_exec(cmd = blast_cmd)
#   return(blast_res)
# }
#
#
# #1c - get_fasta_header ----
#
# ###          function to get fasta names from db based on subjectIDs         ###
# get_fasta_header <- function(id, db_path = "/data/databases/nt/nt") {
#   # id <- "JQ365494.1"
#   command <- "blastdbcmd -db {db_path} -entry {id} -outfmt %t"
#   result <- shell_exec(cmd = command)
#   return(result$stdout)
# }
#
# #1d - get_blast_results ----
# ###          function to get fasta names from db based on subjectIDs         ###
#
# get_blast_results <- function(asv, num_threads = 40) {
#   blast_res <- run_blast(asv, num_threads = num_threads)
#   if (blast_res$status != 0) {
#     rlang::abort(message = "Blast has not run correctly.")
#   }
#   `%>%` <- dplyr::`%>%`
#
#   if (blast_res$stdout == "") {
#     # rlang::inform(glue::glue("Sequence {asv} not found."))
#     df_to_return <- tibble::tibble(`OTU` = asv)
#     return(df_to_return)
#   }
#   blast_table <- blast_res$stdout %>%
#     readr::read_delim(delim = "\t",
#                       col_names = c("query","subject","indentity","length","mismatches","gaps",
#                                     "query start","query end","subject start","subject end",
#                                     "e-value","bitscore","qcovhsp"),
#                       trim_ws = TRUE, comment = "#"
#     )
#
#   blast_table$`subject header` <- purrr::map_chr(blast_table$subject, get_fasta_header)
#   blast_table <- dplyr::relocate(blast_table, `subject header`)
#   blast_table <- tibble::rowid_to_column(blast_table, var = "res")
#
#   # blast_table <- tidyr::pivot_wider(blast_table, id_cols = subject,  names_from = res, values_from = seq_len(ncol(blast_table)), names_glue = "{res}_{.value}")
#   blast_table <- tidyr::pivot_wider(blast_table, names_from = res, values_from = seq_len(ncol(blast_table)), names_glue = "{res}_{.value}")
#
#   blast_table <- blast_table %>%
#     dplyr::mutate(`OTU` = asv) %>%
#     dplyr::relocate(starts_with("3_")) %>%
#     dplyr::relocate(starts_with("2_")) %>%
#     dplyr::relocate(starts_with("1_")) %>%
#     dplyr::relocate(`OTU`)
#   # %>%                                # limpar colunas redundantes. Não está ativado pois estava dando pau
#   #   select(-c("1_res","1_query",
#   #             "2_res","2_query",
#   #             "3_res","3_query"))
#
#   return(blast_table)
# }
#
# #2 - Running BLASTr function with supplied ASVs ----
# #2a - Set parallelization ----
#
# cores_to_be_used <- future::availableCores() - 2 # Usar todos os cores -2 = 78
#
# cat(paste0("This analysis will be performed using \n\t",cores_to_be_used," cores"))
#
# future::plan(future::multisession(workers = cores_to_be_used))
#
# #2b - Run BLAST ----
# blast_res <- furrr::future_map_dfr(asvs_blast, get_blast_results, num_threads = 1, .options = furrr::furrr_options(seed = NULL))
# blast_res <- furrr::future_map_dfr("CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC", get_blast_results, num_threads = 1, .options = furrr::furrr_options(seed = NULL))
#
# #3 - Results ----
# #3a - save blast res to file ----
# readr::write_csv(blast_res, paste0(results_path,"/",prjct_rad,"-asv_blastn_res_",Sys.Date(),".csv"),append = FALSE)
#
#
# #3b - return tibble ----
#
# blast_res  #results table with ASVs identified (or not)
#
#
#
#
