test_that("Check if CLI tool is available", {
  expect_equal(
    check_bin("R"),
    Sys.which("R")
  )
})
