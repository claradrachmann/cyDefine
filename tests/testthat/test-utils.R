test_that("Negated %in% works", {
  expect_true(1 %!in% 2:5)
})

test_that("baseR is installed", {
  expect_true(check_package("base"))
})
