test_that("quaqc errors when bin is NULL", {
  expect_error(
    quaqc(bam.files = "x.bam", bin = NULL),
    "Please specify the location of the quaqc binary"
  )
})

test_that("quaqc errors when peaks is not a character or GRanges", {
  expect_error(
    quaqc(bam.files = "x.bam", peaks = 42L, bin = "quaqc"),
    "Expected a filename or GRanges object for peaks"
  )
})

test_that("quaqc errors when tss is not a character or GRanges", {
  expect_error(
    quaqc(bam.files = "x.bam", tss = 42L, bin = "quaqc"),
    "Expected a filename or GRanges object for tss"
  )
})

test_that("quaqc errors when target.list is not a character or GRanges", {
  expect_error(
    quaqc(bam.files = "x.bam", target.list = 42L, bin = "quaqc"),
    "Expected a filename or GRanges object for target.list"
  )
})

test_that("quaqc errors when blacklist is not a character or GRanges", {
  expect_error(
    quaqc(bam.files = "x.bam", blacklist = 42L, bin = "quaqc"),
    "Expected a filename or GRanges object for blacklist"
  )
})

test_that("quaqc without bam.files returns a character vector of help text", {
  skip_if(Sys.which(getOption("quaqc.bin", "quaqc")) == "", "quaqc binary not on PATH")
  result <- quaqc()
  expect_type(result, "character")
  expect_gt(length(result), 0)
})
