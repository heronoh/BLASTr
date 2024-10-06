test_that("Execute shell commands", {
  testthat::skip_on_os("windows")
  testthat::skip_on_ci()
  expect_equal(
    shell_exec("R -s -q -e 'base::cat(base::version$version.string)'")$stdout,
    base::version$version.string
  )
})
