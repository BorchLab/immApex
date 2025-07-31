# test script for calculateProperty.R - testcases are NOT comprehensive!


test_that("baseline output (mean, Atchley) has correct dimensions & names", {
  seqs <- c("ACDE", "ACDF")                      # L = 4
  
  res  <- calculateProperty(seqs,
                            property.set = "atchleyFactors",
                            summary.fun  = "mean")
  
  expect_type(res, "double")
  expect_equal(dim(res), c(5, 4))                # 5 factors × 4 positions
  expect_equal(rownames(res), paste0("AF", 1:5))
  expect_equal(colnames(res), paste0("Pos.", 1:4))
})

test_that("first column equals simple mean of Atchley values", {
  seqs <- c("AC", "AC")                          # all residues identical
  res  <- calculateProperty(seqs, "atchleyFactors")
  
  Avals <- .builtin_scales$atchleyFactors[ , "A"]        # AF1..AF5 for Alanine
  Cvals <- .builtin_scales$atchleyFactors[ , "C"]        # for checking col-2
  
  expect_equal(res[ , "Pos.1"], Avals, tolerance = 1e-12)
  expect_equal(res[ , "Pos.2"], Cvals, tolerance = 1e-12)
})

test_that("summary.fun = 'sum' is nSeq × 'mean'", {
  seqs <- c("ACDE", "ACDF", "ACDG")
  m    <- calculateProperty(seqs, summary.fun = "mean")
  s    <- calculateProperty(seqs, summary.fun = "sum")
  
  expect_equal(s, m * length(seqs), tolerance = 1e-12)
})

test_that("custom summary function works (max)", {
  seqs <- c("ACDE", "ACDF")
  max_fun <- function(x) max(x, na.rm = TRUE)
  
  res_max <- calculateProperty(seqs, summary.fun = max_fun)
  
  # should be element-wise ≥ mean and identical where only one residue
  res_mean <- calculateProperty(seqs, summary.fun = "mean")
  expect_true(all(res_max >= res_mean))
  expect_equal(res_max[ , "Pos.3"], res_mean[ , "Pos.3"])  # single residue (D)
})

test_that("all transform options behave as expected", {
  seqs <- c("ACDE", "ACDE")
  base <- calculateProperty(seqs, transform = "none")
  
  
  mm <- calculateProperty(seqs, transform = "minmax")
  rng <- apply(mm, 1, range)
  expect_true(all(abs(rng[1, ]) < 1e-12))        # mins ~ 0
  expect_true(all(abs(rng[2, ] - 1) < 1e-12))    # maxs ~ 1
})

test_that("tidy = TRUE returns correct long data.frame", {
  seqs <- c("AC")
  df   <- calculateProperty(seqs, tidy = TRUE)
  
  expect_s3_class(df, "data.frame")
  expect_named(df, c("scale", "position", "value"))
  expect_equal(nrow(df), 5 * 2)                  # 5 scales × L = 2
  expect_equal(df$position, rep(1:2, each = 5))
})

test_that("custom property matrix input works", {
  M <- matrix(runif(40), nrow = 2,
              dimnames = list(paste0("S", 1:2), amino.acids))
  seqs <- c("AC", "AD")
  
  res <- calculateProperty(seqs, property.set = M)
  expect_equal(dim(res), c(2, 2))
  expect_equal(rownames(res), c("S1", "S2"))
})

#  ── Error handling ───────────────────────────────────────────────

test_that("invalid inputs trigger errors", {
  seqs <- c("ACD")
  
  # non-character sequences
  expect_error(calculateProperty(1:3))
  
  # padding symbol duplicates an amino acid
  expect_error(calculateProperty(seqs, padding.symbol = "A"))
  
  # unknown summary.fun keyword
  expect_error(calculateProperty(seqs, summary.fun = "bogus"))
  
  # unknown property set
  expect_error(calculateProperty(seqs, property.set = "NotASet"))
})
