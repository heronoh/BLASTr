testthat::test_that("parse_fasta() reads the sequences and names a missing file", {
  f <- fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta")
  testthat::expect_true(length(parse_fasta(f)) > 0L)
  missing <- fs::file_temp(ext = "fasta")
  testthat::expect_error(parse_fasta(missing), regexp = basename(missing), class = "blastr_fasta_file_not_readable")
})
