test_that("melt_reports errors on non-quaqc input", {
  expect_error(melt_reports(list()), "'report' must be a 'quaqc' class object")
  expect_error(melt_reports("not a report"), "'report' must be a 'quaqc' class object")
})

test_that("melt_reports errors on invalid section argument", {
  expect_error(melt_reports(FIXTURE, "no_such_section"))
})

test_that("bam_stats section returns expected columns", {
  result <- melt_reports(FIXTURE, "bam_stats")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Sample", "Reads", "Count"))
  expect_gt(nrow(result), 0)
})

test_that("overview_unfilt section returns expected columns", {
  result <- melt_reports(FIXTURE, "overview_unfilt")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Sample", "Reads", "Nuclear", "Mitochondrial", "Plastidic"))
  expect_gt(nrow(result), 0)
})

test_that("overview_filt section returns expected columns", {
  result <- melt_reports(FIXTURE, "overview_filt")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Sample", "Reads", "Nuclear", "Mitochondrial", "Plastidic"))
  expect_gt(nrow(result), 0)
})

test_that("nucl_stats section returns expected columns", {
  result <- melt_reports(FIXTURE, "nucl_stats")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Sample", "Reads", "AlignmentSize", "FragmentSize", "MAPQ",
    "ReadDepth", "GCPercent"))
  expect_gt(nrow(result), 0)
})

test_that("nucl_addn section returns expected columns", {
  result <- melt_reports(FIXTURE, "nucl_addn")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Sample", "GenomeCoverage", "AlnNoMAPQ"))
  expect_equal(nrow(result), 1)
})

test_that("peak_stats section returns expected columns", {
  result <- melt_reports(FIXTURE, "peak_stats")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Sample", "PeakCount", "PeakGenomeCov", "FRIP"))
})

test_that("tss_stats section returns expected columns", {
  result <- melt_reports(FIXTURE, "tss_stats")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Sample", "TSSCount", "TES"))
})

test_that("tss_pileup section returns expected columns", {
  result <- melt_reports(FIXTURE, "tss_pileup")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Sample", "Coordinate", "Depth"))
})

test_that("aln_hist section returns expected columns", {
  result <- melt_reports(FIXTURE, "aln_hist")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Sample", "AlnSize", "Count"))
})

test_that("frag_hist section returns expected columns", {
  result <- melt_reports(FIXTURE, "frag_hist")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Sample", "FragSize", "Count"))
})

test_that("gc_hist section returns expected columns", {
  result <- melt_reports(FIXTURE, "gc_hist")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Sample", "GCPercent", "Count"))
})

test_that("depth_hist section returns expected columns", {
  result <- melt_reports(FIXTURE, "depth_hist")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Sample", "ReadDepth", "Count"))
})

test_that("genome section returns expected columns", {
  result <- melt_reports(FIXTURE, "genome")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Sample", "Sequence", "Count", "Size"))
  expect_equal(nrow(result), 5)
})

test_that("use.basename strips directory components from Sample column", {
  with_dir <- melt_reports(FIXTURE, "bam_stats", use.basename = FALSE)
  with_base <- melt_reports(FIXTURE, "bam_stats", use.basename = TRUE)
  expect_equal(with_base$Sample, basename(with_dir$Sample))
})

test_that("normalize.hist = 'proportion' makes counts sum to 1", {
  result <- melt_reports(FIXTURE, "aln_hist", normalize.hist = "proportion")
  expect_equal(sum(result$Count), 1, tolerance = 1e-10)
})

test_that("normalize.hist = 'max' makes the maximum count equal to 1", {
  result <- melt_reports(FIXTURE, "aln_hist", normalize.hist = "max")
  expect_equal(max(result$Count), 1)
})

test_that("normalize.tss options all return numeric Depth values", {
  for (norm in c("no", "bkg", "rpm")) {
    result <- melt_reports(FIXTURE, "tss_pileup", normalize.tss = norm)
    expect_true(is.numeric(result$Depth),
      info = paste("normalize.tss =", norm))
  }
})
