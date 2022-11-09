test_that("Execute shell commands", {
  expect_equal(
    shell_exec("R -s -q -e 'base::cat(base::version$version.string)'")$stdout,
    base::version$version.string
  )
})
