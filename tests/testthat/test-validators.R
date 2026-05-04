test_that("validate_quaqc returns TRUE for a valid quaqc object", {
  expect_true(quaqcr:::validate_quaqc(FIXTURE))
})

test_that("validate_quaqc errors on non-quaqc input", {
  expect_error(quaqcr:::validate_quaqc(list()), "Object is not of class 'quaqc'")
  expect_error(quaqcr:::validate_quaqc(42L), "Object is not of class 'quaqc'")
})

test_that("validate_quaqc_report returns TRUE for a valid quaqc_report object", {
  expect_true(quaqcr:::validate_quaqc_report(FIXTURE$reports[[1]]))
})

test_that("validate_quaqc_report errors on non-quaqc_report input", {
  expect_error(
    quaqcr:::validate_quaqc_report(list()),
    "Object is not of class 'quaqc_report'"
  )
  expect_error(
    quaqcr:::validate_quaqc_report(FIXTURE),
    "Object is not of class 'quaqc_report'"
  )
})

test_that("identical_params returns TRUE when comparing a report with itself", {
  rep <- FIXTURE$reports[[1]]
  expect_true(quaqcr:::identical_params(rep, rep))
})

test_that("identical_params errors on non-quaqc_report inputs", {
  expect_error(
    quaqcr:::identical_params(list(), list()),
    "'x' and 'y' must be 'quaqc_report' objects"
  )
  expect_error(
    quaqcr:::identical_params(FIXTURE$reports[[1]], list()),
    "'x' and 'y' must be 'quaqc_report' objects"
  )
})

test_that("remove_fails keeps successful reports unchanged", {
  result <- quaqcr:::remove_fails(FIXTURE)
  expect_equal(length(result$reports), 1)
  expect_equal(result$metadata$samples, 1)
})

test_that("remove_fails removes failed reports and updates sample count", {
  failed <- FIXTURE
  failed$reports[[1]]$success <- FALSE
  result <- quaqcr:::remove_fails(failed)
  expect_equal(length(result$reports), 0)
  expect_equal(result$metadata$samples, 0)
})
