test_that("footprint errors when json = NULL", {
  expect_error(
    footprint("p.bed", "x.bam", json = NULL),
    "'json' cannot be NULL"
  )
})

test_that("pileup errors when json = NULL", {
  expect_error(
    pileup("p.bed", "x.bam", json = NULL),
    "'json' cannot be NULL"
  )
})

test_that("pileup warns when tss is provided", {
  # pileup() has no tss.* formals, so tss goes cleanly into ... and the
  # warning fires before the subsequent quaqc() call errors out.
  expect_warning(
    try(pileup("p.bed", "x.bam", tss = "t.bed"), silent = TRUE),
    "'tss' cannot be set"
  )
})
