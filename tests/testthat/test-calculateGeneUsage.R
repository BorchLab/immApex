# test script for calculateGeneUsage.R - testcases are NOT comprehensive!

df_chr <- data.frame(
  V = c("V1", "V1", "V2", "V3"),
  J = c("J1", "J2", "J1", "J1"),
  stringsAsFactors = FALSE
)

df_fac <- within(df_chr, {
  V <- factor(V, levels = c("V1","V2","V3"))
  J <- factor(J, levels = c("J1","J2"))
})

# ---------------------------------------------------------------------------
# single-locus summaries 
# ---------------------------------------------------------------------------
test_that("single-locus 'proportion' returns named numeric vector", {
  res <- calculateGeneUsage(df_chr, loci = "V")
  expect_type(res, "double")
  expect_named(res, c("V1","V2","V3"))
  expect_equal(sum(res), 1, tolerance = 1e-12)
  expect_equal(as.numeric(res["V1"]), 0.5)
})

test_that("single-locus 'count' and 'percent' scale correctly", {
  cnt <- calculateGeneUsage(df_chr, loci = "V", summary = "count")
  pct <- calculateGeneUsage(df_chr, loci = "V", summary = "percent")
  expect_equal(sum(cnt), nrow(df_chr))
  expect_equal(sum(pct), 100, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# paired-locus summaries 
# ---------------------------------------------------------------------------
test_that("paired-locus 'proportion' returns matrix with correct dims", {
  prop <- calculateGeneUsage(df_chr, loci = c("V","J"))
  expect_true(is.matrix(prop))
  expect_equal(dim(prop), c(length(unique(df_chr$V)), length(unique(df_chr$J))))
  expect_equal(sum(prop), 1, tolerance = 1e-12)
  # V1-J1 cell should be 1/4
  expect_equal(prop["V1","J1"], 0.25)
})

test_that("levels argument pads missing combinations with zeros", {
  levs <- list(c("V1","V2","V3","V4"), c("J1","J2","J3"))
  mat  <- calculateGeneUsage(df_chr, loci = c("V","J"), levels = levs,
                             summary = "count")
  expect_equal(dim(mat), c(4L, 3L))
  expect_equal(mat["V4","J3"], 0)
  # row/column totals still equal overall n
  expect_equal(sum(mat), nrow(df_chr))
})

# ---------------------------------------------------------------------------
# factor input behaves like character 
# ---------------------------------------------------------------------------

test_that("factor columns are handled correctly", {
  res1 <- calculateGeneUsage(df_chr, loci = "V")
  res2 <- calculateGeneUsage(df_fac, loci = "V")
  expect_equal(res1, res2, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# error handling 
# ---------------------------------------------------------------------------
test_that("invalid inputs raise informative errors", {
  # loci not present
  expect_error(calculateGeneUsage(df_chr, loci = "X"))
  # more than two loci
  expect_error(calculateGeneUsage(df_chr, loci = c("V","J","X")))
  # non-data.frame input
  expect_error(calculateGeneUsage(list(V="V1"), loci = "V"))
  # bad summary argument
  expect_error(calculateGeneUsage(df_chr, loci = "V", summary = "bogus"))
  # mismatched levels length (paired but only one levels vector supplied)
  bad_levels <- list(c("V1","V2"))
  expect_error(calculateGeneUsage(df_chr, loci = c("V","J"), levels = bad_levels))
})
