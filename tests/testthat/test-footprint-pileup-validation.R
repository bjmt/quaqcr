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

test_that("pileup errors when tss is provided", {
  expect_error(
    pileup("p.bed", "x.bam", tss = "t.bed"),
    "'tss' cannot be set"
  )
})
