test_that("parse_quaqc_file returns a quaqc object with the right structure", {
  expect_s3_class(FIXTURE, "quaqc")
  expect_named(FIXTURE, c("metadata", "reports"))
})

test_that("quaqc metadata has expected fields and values", {
  meta <- FIXTURE$metadata
  expect_named(meta, c("version", "title", "date", "args", "samples"))
  expect_equal(meta$samples, 1L)
  expect_type(meta$version, "character")
  expect_type(meta$date, "character")
})

test_that("quaqc reports list has one quaqc_report", {
  expect_length(FIXTURE$reports, 1)
  expect_s3_class(FIXTURE$reports[[1]], "quaqc_report")
})

test_that("the fixture quaqc_report is successful and has all sub-structures", {
  rep <- FIXTURE$reports[[1]]
  expect_true(rep$success)
  expect_type(rep$sample, "character")
  expect_true(is.matrix(rep$genome))
  expect_true(is.matrix(rep$unfiltered$effective))
  expect_s3_class(rep$filtered$nuclear$stats, "data.frame")
  expect_named(rep$params, c("integer", "boolean"))
})

test_that("parse_quaqc errors on a non-list input", {
  expect_error(parse_quaqc("not a list"))
  expect_error(parse_quaqc(42L))
})

test_that("parse_quaqc errors when JSON names are not the expected keys", {
  bad_json <- list(bad_key = 1, another_bad = 2)
  expect_error(parse_quaqc(bad_json))
})

test_that("parse_quaqc_file errors on a non-length-1 non-character vector", {
  expect_error(
    parse_quaqc_file(c(1L, 2L)),
    "'json.file' must be a length 1 character vector"
  )
})

test_that("parse_quaqc_file round-trips through parse_quaqc", {
  fixture.path <- system.file("extdata", "report.json.gz", package = "quaqcr")
  f <- gzfile(fixture.path, "rt")
  json <- jsonlite::fromJSON(readLines(f), simplifyDataFrame = FALSE)
  close(f)
  result <- parse_quaqc(json)
  expect_s3_class(result, "quaqc")
  expect_equal(result$metadata$samples, FIXTURE$metadata$samples)
})
