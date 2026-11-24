testthat::test_that("parse_blast_version reads the version of `blastn -version`", {
  testthat::expect_equal(parse_blast_version("blastn: 2.17.0+\n Package: blast 2.17.0, build Aug 22 2026"), "2.17.0")
  testthat::expect_equal(parse_blast_version("blastn: 2.12.0+"), "2.12.0")
  testthat::expect_true(is.na(parse_blast_version("")))
})

# A fake BLAST+ program that answers `-version` with `version`.
local_fake_blast <- function(version, env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  bin <- file.path(dir, "blastn")
  writeLines(c("#!/bin/sh", sprintf("echo 'blastn: %s+'", version)), bin)
  Sys.chmod(bin, "0755")
  withr::local_envvar(PATH = dir, .local_envir = env)
  dir
}

testthat::test_that("without the conda environment, BLAST+ of the system is used only if recent", {
  testthat::skip_on_os("windows")
  local_fake_blast("2.12.0")
  testthat::expect_error(resolve_blast_bin("blastn", "blastr-no-such-env"), class = "blastr_error_blast_not_found")
  testthat::expect_error(resolve_blast_bin("blastn", "blastr-no-such-env"), regexp = "2.12.0")
  testthat::expect_false(has_blast("blastr-no-such-env"))
})

testthat::test_that("BLAST+ of the system, recent enough, is used when the environment has none", {
  testthat::skip_on_os("windows")
  dir <- local_fake_blast("2.16.0")
  b <- resolve_blast_bin("blastn", "blastr-no-such-env")
  testthat::expect_equal(b$source, "system")
  testthat::expect_equal(b$version, "2.16.0")
  testthat::expect_equal(b$path, normalizePath(file.path(dir, "blastn")))
  testthat::expect_true(has_blast("blastr-no-such-env"))
  testthat::expect_message(check_cmd("blastn", env_name = "blastr-no-such-env"), "of the system")
  testthat::expect_true(is.na(env_blast_version("blastr-no-such-env")))
})

testthat::test_that("nothing found: an error that says to run install_dependencies()", {
  withr::local_envvar(PATH = withr::local_tempdir())
  testthat::expect_error(blast_bin_path("blastn", "blastr-no-such-env"), regexp = "install_dependencies")
})

testthat::test_that("the conda environment is used first", {
  skip_if_no_blast()
  testthat::skip_if(is.na(env_blast_version("blastr-blast-env")), "no conda environment")
  testthat::expect_equal(resolve_blast_bin("blastn")$source, "conda")
})
