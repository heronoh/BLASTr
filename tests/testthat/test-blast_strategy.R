fake_db <- function(dir, name, sizes) {
  fs::dir_create(dir)
  for (i in seq_along(sizes)) {
    f <- fs::path(dir, sprintf("%s.%02d.nsq", name, i - 1L))
    con <- file(f, "wb")
    seek(con, sizes[[i]] - 1)
    writeBin(as.raw(0), con)
    close(con)
  }
  file.create(fs::path(dir, paste0(name, ".nal")))
  as.character(fs::path(dir, name))
}

testthat::test_that("small database: one single-threaded process per core", {
  d <- withr::local_tempdir()
  db <- fake_db(d, "small", c(1000, 1000))
  s <- blast_strategy(n_seqs = 100, db_path = db, total_cores = 8, ram_gb = 16)
  testthat::expect_equal(s$threads, 1L)
  testthat::expect_equal(s$processes, 8L)
  testthat::expect_equal(s$n_batches, 8L)
  testthat::expect_equal(s$db_volumes, 2L)
})

testthat::test_that("large database: threads split the database", {
  d <- withr::local_tempdir()
  db <- fake_db(d, "big", rep(1e9, 4))
  s <- blast_strategy(n_seqs = 100, db_path = db, total_cores = 8, ram_gb = 64, small_db_gb = 1)
  testthat::expect_equal(s$mt_mode, 2L)
  testthat::expect_equal(s$threads, 4L) # limited by the 4 volumes
  testthat::expect_equal(s$processes, 2L)
  testthat::expect_true(s$fits_ram)

  s <- blast_strategy(n_seqs = 100, db_path = db, total_cores = 8, ram_gb = 2)
  testthat::expect_false(s$fits_ram)
  testthat::expect_equal(s$processes, 1L)
  testthat::expect_equal(s$threads, 8L)
  # one batch: each batch would read the whole database from disk again
  s <- blast_strategy(n_seqs = 3000, db_path = db, total_cores = 8, ram_gb = 2)
  testthat::expect_equal(s$n_batches, 1L)
  s <- blast_strategy(n_seqs = 3000, db_path = db, total_cores = 8, ram_gb = 2, batch_size = 1000)
  testthat::expect_equal(s$n_batches, 3L)
})

testthat::test_that("several databases: sizes added, smallest first", {
  d <- withr::local_tempdir()
  big <- fake_db(d, "big", rep(1e9, 2))
  small <- fake_db(d, "small", 1000)
  s <- blast_strategy(n_seqs = 10, db_path = c(big, small), total_cores = 2, ram_gb = 64)
  testthat::expect_equal(s$db_path, c(small, big))
  testthat::expect_equal(s$db_gb, 2 + 1e-6, tolerance = 1e-6)
  s2 <- blast_strategy(n_seqs = 10, db_path = paste(big, small), total_cores = 2, ram_gb = 64)
  testthat::expect_equal(s2$db_path, c(small, big))
})

testthat::test_that("manual values override the automatic strategy", {
  d <- withr::local_tempdir()
  db <- fake_db(d, "small", 1000)
  s <- blast_strategy(n_seqs = 100, db_path = db, total_cores = 8, ram_gb = 16,
                      num_threads = 4, batch_size = 10, mt_mode = 1)
  testthat::expect_equal(s$threads, 4L)
  testthat::expect_equal(s$processes, 2L)
  testthat::expect_equal(s$batch_size, 10L)
  testthat::expect_equal(s$n_batches, 10L)
  testthat::expect_equal(s$mt_mode, 1L)
})

testthat::test_that("batches are capped and never fewer than processes", {
  d <- withr::local_tempdir()
  db <- fake_db(d, "small", 1000)
  s <- blast_strategy(n_seqs = 3000, db_path = db, total_cores = 4, ram_gb = 16, max_batch_size = 500)
  testthat::expect_lte(s$batch_size, 500L)
  testthat::expect_equal(s$n_batches %% 4L, 0L)
  s <- blast_strategy(n_seqs = 2, db_path = db, total_cores = 8, ram_gb = 16)
  testthat::expect_equal(s$processes, 2L)
})

testthat::test_that("core_nt-like database in RAM: one process, few batches", {
  # 92 volumes of ~3.2 GB, as core_nt in the benchmark (sparse files)
  d <- withr::local_tempdir()
  db <- fake_db(d, "core", rep(3.2e9, 92))
  s <- blast_strategy(n_seqs = 2000, db_path = db, total_cores = 32, ram_gb = 500)
  testthat::expect_equal(s$processes, 1L)
  testthat::expect_equal(s$threads, 32L)
  testthat::expect_equal(s$mt_mode, 2L)
  testthat::expect_equal(s$n_batches, 2L)
  s <- blast_strategy(n_seqs = 100, db_path = db, total_cores = 32, ram_gb = 500)
  testthat::expect_equal(s$n_batches, 1L)
  # more than 32 cores: more processes of 32 threads, not more threads
  s <- blast_strategy(n_seqs = 2000, db_path = db, total_cores = 64, ram_gb = 500)
  testthat::expect_equal(s$processes, 2L)
  testthat::expect_equal(s$threads, 32L)
  testthat::expect_equal(s$n_batches, 2L)
  s <- blast_strategy(n_seqs = 2000, db_path = db, total_cores = 64, ram_gb = 500, max_threads = 64)
  testthat::expect_equal(s$processes, 1L)
})
