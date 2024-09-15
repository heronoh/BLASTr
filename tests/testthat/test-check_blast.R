# test_that("Check if BLAST CLI is available", {
#  expect_equal(class(check_bin("blastn")), "character")
# })

test_that("Error if BLAST not on PATH", {
  withr::local_envvar(list(`PATH` = ""))
  expect_error(
    check_bin("blastn"),
    class = "blastr_bin_not_found"
  )
})
