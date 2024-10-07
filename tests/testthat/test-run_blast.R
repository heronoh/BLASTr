test_that("parallel_blast works", {
  ASVs_test <- readLines(fs::path_package("BLASTr", "extdata", "asvs_test", ext = "txt"))
  db_path <- fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta")
  blast_res <- parallel_blast(
    asvs = ASVs_test,
    db_path = db_path,
    out_file = NA,
    out_RDS = NA,
    total_cores = 4,
    perc_id = 80,
    num_threads = 1,
    perc_qcov_hsp = 80,
    num_alignments = 4,
    blast_type = "blastn",
    verbose = FALSE
  )

  testthat::expect_s3_class(blast_res, "tbl_df")
})


test_that("parallel_blast works", {
  ASVs_test <- readLines(fs::path_package("BLASTr", "extdata", "asvs_test", ext = "txt"))
  db_path <- fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta")

  blast_res <- run_blast(
    asv = ASVs_test[1],
    db_path = db_path
  )

  testthat::expect_equal(blast_res$status, 0)
})
