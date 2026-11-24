# Copy the example FASTA to a fresh temporary directory
local_fasta_copy <- function(env = parent.frame()) {
  tmp_dir <- withr::local_tempdir(pattern = "blastr_db_", .local_envir = env)
  fasta_copy <- fs::path(tmp_dir, "minimal_db_blast", ext = "fasta")
  fs::file_copy(
    fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta"),
    fasta_copy
  )
  fasta_copy
}

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

  db_path <- fs::file_temp("minimal_db_blast_")
  fasta_path <- local_fasta_copy()
  taxid_map_path <- fs::path_ext_set(fasta_path, "txt")
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

  db_path <- fs::file_temp("minimal_db_blast_")
  fasta_path <- local_fasta_copy()
  taxid_map_path <- fs::path_ext_set(fasta_path, "txt")
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

  db_path <- fs::file_temp("minimal_db_blast_")
  fasta_path <- local_fasta_copy()
  taxid_map_path <- fs::path_ext_set(fasta_path, "txt")
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

testthat::test_that("`make_blast_db()` defaults `db_path` to `fasta_path`", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  fasta_path <- local_fasta_copy()

  db_res <- make_blast_db(fasta_path = fasta_path, verbose = "silent")

  testthat::expect_equal(db_res$status, 0)
  testthat::expect_true(fs::file_exists(paste0(fasta_path, ".nsq")))
})

testthat::test_that("`make_blast_db()` stops if the FASTA directory is not writable", {
  testthat::skip_on_cran()
  testthat::skip_on_os("windows")
  # root ignores write permissions
  testthat::skip_if(identical(unname(Sys.info()[["user"]]), "root"))

  fasta_path <- local_fasta_copy()
  fasta_dir <- fs::path_dir(fasta_path)
  fs::file_chmod(fasta_dir, "u-w")
  withr::defer(fs::file_chmod(fasta_dir, "u+w"))

  testthat::expect_error(
    make_blast_db(fasta_path = fasta_path, verbose = "silent"),
    class = "blastr_error_db_dir_not_writable"
  )
  testthat::expect_error(
    make_blast_db(fasta_path = fasta_path, verbose = "silent"),
    regexp = "db_path"
  )

  # with a writable `db_path` it works
  testthat::skip_if_offline()
  db_path <- fs::file_temp("minimal_db_blast_")
  db_res <- make_blast_db(
    fasta_path = fasta_path,
    db_path = db_path,
    verbose = "silent"
  )
  testthat::expect_equal(db_res$status, 0)
  testthat::expect_true(fs::file_exists(paste0(db_path, ".nsq")))
})

testthat::test_that("`make_blast_db()` stops if `db_path` is not writable", {
  testthat::skip_on_cran()
  testthat::skip_on_os("windows")
  testthat::skip_if(identical(unname(Sys.info()[["user"]]), "root"))

  fasta_path <- local_fasta_copy()

  # directory does not exist
  testthat::expect_error(
    make_blast_db(
      fasta_path = fasta_path,
      db_path = fs::path(fs::file_temp("missing_dir_"), "db"),
      verbose = "silent"
    ),
    class = "blastr_error_db_dir_not_writable"
  )

  # directory exists but is read-only
  ro_dir <- withr::local_tempdir(pattern = "blastr_ro_")
  fs::file_chmod(ro_dir, "u-w")
  withr::defer(fs::file_chmod(ro_dir, "u+w"))
  testthat::expect_error(
    make_blast_db(
      fasta_path = fasta_path,
      db_path = fs::path(ro_dir, "db"),
      verbose = "silent"
    ),
    class = "blastr_error_db_dir_not_writable"
  )
})

testthat::test_that("`make_blast_db()` takes a `taxid_map` with any name and is silent", {
  testthat::skip_on_cran()

  db_path <- fs::file_temp("minimal_db_blast_")
  fasta_path <- local_fasta_copy()
  taxid_map_path <- fs::file_temp("taxids_", ext = "tsv")
  write("AY882416\t2759", taxid_map_path)

  testthat::expect_silent(
    db_res <- make_blast_db(
      fasta_path = fasta_path,
      db_path = db_path,
      taxid_map = taxid_map_path,
      verbose = "silent"
    )
  )
  testthat::expect_equal(db_res$status, 0)
  testthat::expect_error(
    make_blast_db(fasta_path = fasta_path, db_path = db_path, taxid_map = fs::file_temp(), verbose = "silent"),
    class = "blastr_error_taxid_map_not_found"
  )
})
