test_that("print.quaqc output matches snapshot", {
  expect_snapshot(print(FIXTURE))
})

test_that("print.quaqc_report output matches snapshot", {
  expect_snapshot(print(FIXTURE$reports[[1]]))
})
