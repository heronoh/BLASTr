testthat::test_that("`make_blast_db()` works", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()
  fasta_path <- fs::path_package(
    "BLASTr",
    "extdata",
    "minimal_db_blast",
    ext = "fasta"
  )
  db_path <- fs::file_temp("minimal_db_blast_")

  db_res <- make_blast_db(
    fasta_path = fasta_path,
    db_path = db_path,
    db_type = "nucl",
    verbose = "silent"
  )

  testthat::expect_equal(db_res$status, 0)
  testthat::expect_true(fs::file_exists(paste0(db_path, ".nsq")))
})

testthat::test_that("`make_blast_db()` fails with wrong db_type", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  fasta_path <- fs::path_package(
    "BLASTr",
    "extdata",
    "minimal_db_blast",
    ext = "fasta"
  )
  db_path <- fs::file_temp("minimal_db_blast_")

  testthat::expect_error(
    make_blast_db(
      fasta_path = fasta_path,
      db_path = db_path,
      db_type = "INVALID_DB_TYPE",
      verbose = "silent"
    ),
    class = "blastr_error_make_blast_db"
  )
})

testthat::test_that("`make_blast_db()` works with `taxid_map`", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  fasta_path <- fs::path_package(
    "BLASTr",
    "extdata",
    "minimal_db_blast",
    ext = "fasta"
  )
  db_path <- fs::file_temp("minimal_db_blast_")
  taxid_map_path <- fs::file_temp("taxid_map_", ext = "tsv")
  write("AY882416	2759", taxid_map_path)

  db_res <- make_blast_db(
    fasta_path = fasta_path,
    db_path = db_path,
    db_type = "nucl",
    taxid_map = taxid_map_path,
    verbose = "silent"
  )

  testthat::expect_equal(db_res$status, 0)
  testthat::expect_true(fs::file_exists(paste0(db_path, ".nsq")))
})

testthat::test_that("`make_blast_db()` works with `parse_seqids`", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  fasta_path <- fs::path_package(
    "BLASTr",
    "extdata",
    "minimal_db_blast",
    ext = "fasta"
  )
  db_path <- fs::file_temp("minimal_db_blast_")
  taxid_map_path <- fs::file_temp("taxid_map_", ext = "tsv")
  write("AY882416	2759", taxid_map_path)

  db_res <- make_blast_db(
    fasta_path = fasta_path,
    db_path = db_path,
    db_type = "nucl",
    taxid_map = taxid_map_path,
    verbose = "silent"
  )

  testthat::expect_equal(db_res$status, 0)
  testthat::expect_true(fs::file_exists(paste0(db_path, ".nsq")))
})

testthat::test_that("`make_blast_db()` output has valid taxid", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()
  # TODO: @luciorq Remove `skip` when example fasta actually
  # + have taxid in the header
  testthat::skip("The example fasta does not have taxid in the header")

  fasta_path <- fs::path_package(
    "BLASTr",
    "extdata",
    "minimal_db_blast",
    ext = "fasta"
  )
  db_path <- fs::file_temp("minimal_db_blast_")
  taxid_map_path <- fs::file_temp("taxid_map_", ext = "tsv")
  base::write("AY882416	2759", taxid_map_path)

  db_res <- make_blast_db(
    fasta_path = fasta_path,
    db_path = db_path,
    db_type = "nucl",
    taxid_map = taxid_map_path,
    verbose = "silent"
  )

  testthat::expect_true(fs::file_exists(fs::path(db_path)))

  testthat::expect_equal(db_res$status, 0)
  testthat::expect_true(fs::file_exists(paste0(db_path, ".nsq")))

  # Test parallel_get_tax() on the created database
  tax_res <- parallel_get_tax(
    db_path = db_path,
    taxids = c(2759)
  )

  testthat::expect_true(is.data.frame(tax_res))
  testthat::expect_true(nrow(tax_res) > 0)
})
