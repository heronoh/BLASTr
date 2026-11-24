testthat::test_that("make_taxonomy_db and get_tax build the lineage", {
  dir <- tiny_taxdump(withr::local_tempdir())
  db <- make_taxonomy_db(file.path(withr::local_tempdir(), "tax.sqlite"), taxdump = dir) |>
    suppressMessages()
  res <- get_tax(c("9606", "9598", "9606"), taxonomy_db = db)
  testthat::expect_equal(res$taxid, c("9606", "9598"))
  testthat::expect_equal(res$sci_name, c("Homo sapiens", "Pan troglodytes"))
  testthat::expect_equal(res$rank, c("species", "species"))
  # superkingdom fills domain (taxdumps before 2025)
  testthat::expect_equal(res$domain, c("Eukaryota", "Eukaryota"))
  testthat::expect_equal(res$genus, c("Homo", "Pan"))
  testthat::expect_true(all(is.na(res$class)))
  testthat::expect_equal(res$lineage[1], "Eukaryota; Metazoa; Chordata; Primates; Hominoidea; Hominidae; Homo; Homo sapiens")
  testthat::expect_equal(res$lineage_taxids[1], "2759; 33208; 7711; 9443; 314295; 9604; 9605; 9606")
})

testthat::test_that("get_tax resolves merged taxids and warns on missing ones", {
  dir <- tiny_taxdump(withr::local_tempdir())
  db <- make_taxonomy_db(file.path(withr::local_tempdir(), "tax.sqlite"), taxdump = dir) |>
    suppressMessages()
  testthat::expect_warning(res <- get_tax(c("9999", "123"), taxonomy_db = db), "not found")
  testthat::expect_equal(res$found, c(TRUE, FALSE))
  testthat::expect_equal(res$merged_into, c("9606", NA))
  testthat::expect_equal(res$species, c("Homo sapiens", NA))
})

testthat::test_that("get_tax takes the taxids of a BLASTr result", {
  dir <- tiny_taxdump(withr::local_tempdir())
  db <- make_taxonomy_db(file.path(withr::local_tempdir(), "tax.sqlite"), taxdump = dir) |>
    suppressMessages()
  withr::local_options(BLASTr.taxonomy_db = db)
  long <- data.frame(Sequence = c("A", "A", "B"), staxid = c("9606", "0", "9598;9606"))
  testthat::expect_equal(get_tax(long)$taxid, c("9606", "9598"))
  wide <- data.frame(Sequence = "A", `1_staxid` = "9598", `2_staxid` = NA, check.names = FALSE)
  testthat::expect_equal(get_tax(wide)$taxid, "9598")
  testthat::expect_equal(nrow(get_tax(character(0L))), 0L)
  testthat::expect_error(get_tax(data.frame(x = 1)), class = "blastr_no_taxids_error")
})

testthat::test_that("get_tax builds a missing database from the taxdump", {
  dir <- tiny_taxdump(withr::local_tempdir())
  db <- file.path(withr::local_tempdir(), "sub", "tax.sqlite")
  res <- get_tax("9606", taxonomy_db = db, taxdump = dir) |> suppressMessages()
  testthat::expect_true(file.exists(db))
  testthat::expect_equal(res$genus, "Homo")
})

testthat::test_that("get_tax keeps the database of a search in the cache, by database date", {
  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir(), BLASTR_TAXONOMY_DB = "")
  withr::local_options(BLASTr.taxonomy_db = NULL)
  dir <- tiny_taxdump(withr::local_tempdir())
  res <- data.frame(Sequence = "A", staxid = "9598")
  attr(res, "BLASTr_metadata") <- list(db_dates = c("2026-07-20T06:00:00", "2026-11-12T01:19:00"))
  tax <- get_tax(res, taxdump = dir) |> suppressMessages()
  testthat::expect_equal(tax$species, "Pan troglodytes")
  testthat::expect_true(file.exists(get_blastr_cache("taxonomy", "taxonomy_2026-11-12.sqlite")))
  # taxids alone: the newest database in the cache, no download
  testthat::expect_equal(get_tax("9606")$genus, "Homo")
})

testthat::test_that("taxdump_url picks the first NCBI archive after the database", {
  testthat::expect_match(taxdump_url("2024-03-15"), "taxdmp_2024-04-01.zip$")
  testthat::expect_match(taxdump_url("2024-03-01"), "taxdmp_2024-03-01.zip$")
  testthat::expect_match(taxdump_url(format(Sys.Date())), "taxdump.tar.gz$")
  testthat::expect_match(taxdump_url(NA), "taxdump.tar.gz$")
})

testthat::test_that("extract_taxids reads `Taxids at best` of wide results (and its old name)", {
  wide <- data.frame(`1_staxid` = "9606", `Taxids at best` = "9606;9598", check.names = FALSE)
  testthat::expect_equal(extract_taxids(wide), c("9606", "9598"))
  names(wide)[2] <- "Taxids at best identity"
  testthat::expect_equal(extract_taxids(wide), c("9606", "9598"))
})

testthat::test_that("get_tax() downloads the taxdump only when asked", {
  db <- fs::path(withr::local_tempdir(), "taxonomy.sqlite")
  # not interactive, no taxdump, no download = TRUE: stops, nothing downloaded
  testthat::expect_error(get_tax("9606", taxonomy_db = db), class = "blastr_error_no_taxonomy_db")
  testthat::expect_error(get_tax("9606", taxonomy_db = db), regexp = "download = TRUE")
  testthat::expect_false(fs::file_exists(db))
  # a taxdump given is a request: the database is built from it
  dir <- tiny_taxdump(withr::local_tempdir())
  res <- suppressMessages(get_tax("9606", taxonomy_db = db, taxdump = dir))
  testthat::expect_equal(res$sci_name, "Homo sapiens")
  testthat::expect_true(fs::file_exists(db))
})

testthat::test_that("has_taxonomy_db() checks the database get_tax() would use, creating nothing", {
  cache_dir <- withr::local_tempdir()
  withr::local_envvar(R_USER_CACHE_DIR = cache_dir, BLASTR_TAXONOMY_DB = "")
  withr::local_options(BLASTr.taxonomy_db = NULL)
  # empty cache: FALSE, no directory created, no download
  x <- has_taxonomy_db()
  testthat::expect_false(x)
  testthat::expect_false(fs::dir_exists(get_blastr_cache()))
  testthat::expect_equal(attr(x, "path"), taxonomy_db_path())

  dir <- tiny_taxdump(withr::local_tempdir())
  db <- make_taxonomy_db(fs::path(withr::local_tempdir(), "tax.sqlite"), taxdump = dir) |> suppressMessages()
  testthat::expect_true(has_taxonomy_db(db))
  testthat::expect_equal(attr(has_taxonomy_db(db), "path"), as.character(db))
  # the option and the environment variable, as in get_tax()
  withr::with_options(list(BLASTr.taxonomy_db = db), testthat::expect_true(has_taxonomy_db()))
  withr::with_envvar(c(BLASTR_TAXONOMY_DB = db), testthat::expect_true(has_taxonomy_db()))

  # a file that is not SQLite, and an SQLite file without the taxonomy tables
  junk <- fs::path(withr::local_tempdir(), "junk.sqlite")
  writeLines("not a database", junk)
  testthat::expect_false(has_taxonomy_db(junk))
  empty <- fs::path(withr::local_tempdir(), "empty.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), empty)
  DBI::dbWriteTable(con, "nodes", data.frame(id = 1L))
  DBI::dbDisconnect(con)
  testthat::expect_false(has_taxonomy_db(empty))
})

testthat::test_that("has_taxonomy_db(db_path =) looks for the database of the BLAST database date", {
  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir(), BLASTR_TAXONOMY_DB = "")
  withr::local_options(BLASTr.taxonomy_db = NULL)
  blast_db <- fs::path(withr::local_tempdir(), "core_nt")
  writeLines('{\n  "last-updated": "2026-11-12T01:19:00"\n}', paste0(blast_db, ".njs"))
  x <- has_taxonomy_db(db_path = blast_db)
  testthat::expect_false(x)
  testthat::expect_equal(fs::path_file(attr(x, "path")), "taxonomy_2026-11-12.sqlite")
  # the database get_tax() builds for the results of a search in it
  res <- data.frame(Sequence = "A", staxid = "9606")
  attr(res, "BLASTr_metadata") <- list(db_dates = db_dates(blast_db))
  get_tax(res, taxdump = tiny_taxdump(withr::local_tempdir())) |> suppressMessages()
  testthat::expect_true(has_taxonomy_db(db_path = blast_db))
  # a newer taxonomy database in the cache does not count for an older BLAST database
  fs::file_move(attr(x, "path"), fs::path(fs::path_dir(attr(x, "path")), "taxonomy_2026-12-01.sqlite"))
  testthat::expect_false(has_taxonomy_db(db_path = blast_db))
  testthat::expect_true(has_taxonomy_db())
})

testthat::test_that("taxonomy_db = NA or \"\" is not set, like NULL", {
  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir(), BLASTR_TAXONOMY_DB = "")
  withr::local_options(BLASTr.taxonomy_db = NULL)
  dir <- tiny_taxdump(withr::local_tempdir())
  db <- get_blastr_cache("taxonomy", "taxonomy_2026-11-12.sqlite")
  fs::dir_create(fs::path_dir(db))
  make_taxonomy_db(db, taxdump = dir) |> suppressMessages()
  for (x in list(NULL, NA, NA_character_, "")) {
    testthat::expect_true(has_taxonomy_db(x))
    testthat::expect_equal(attr(has_taxonomy_db(x), "path"), as.character(db))
    testthat::expect_equal(get_tax("9606", taxonomy_db = x)$genus, "Homo")
  }
  # NA in the option falls back to the environment variable
  withr::local_options(BLASTr.taxonomy_db = NA_character_)
  withr::local_envvar(BLASTR_TAXONOMY_DB = db)
  testthat::expect_equal(attr(has_taxonomy_db(NA), "path"), as.character(db))
})
