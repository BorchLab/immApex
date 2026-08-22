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
  
  # should be element-wise >= mean and identical where only one residue
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

test_that(".aa.property.matrix always returns canonical amino-acid column order", {
  skip_if_not_installed("Peptides")

  # Regression test: sequenceEncoder()/sequenceDecoder() assume column i of
  # this helper's output IS canonical-order amino acid i, positionally (see
  # sequenceEncoder.R's roxygen contract). Peptides::AAdata stores some
  # scales in a different order (MSWHIM/ProtFP: alphabetical; crucianiProperties:
  # E/Q swapped), so this helper must always reorder before returning -- not
  # rely on each caller to do it (calculateProperty() used to be the only one
  # that did).
  peptides_sets <- c("crucianiProperties", "FASGAI", "kideraFactors", "MSWHIM",
                     "ProtFP", "stScales", "tScales", "VHSE", "zScales")
  for (key in peptides_sets) {
    expect_equal(colnames(.aa.property.matrix(key)), amino.acids, info = key)
  }
  expect_equal(colnames(.aa.property.matrix("atchleyFactors")), amino.acids)
})

test_that(".aa.property.matrix aligns columns to a custom sequence.dictionary", {
  skip_if_not_installed("Peptides")

  # The 20 canonical residues in a non-canonical order: alignment must follow
  # the dictionary, not silently fall back to `amino.acids`.
  dict <- rev(amino.acids)
  m <- .aa.property.matrix("MSWHIM", dict)
  expect_identical(colnames(m), dict)
  expect_equal(unname(m[, "R"]), unname(Peptides::mswhimScores("R")[[1]]))

  # A restricted alphabet is a legitimate subset, not an error.
  expect_identical(colnames(.aa.property.matrix("MSWHIM", c("A", "C", "R"))),
                   c("A", "C", "R"))

  # A residue the scale has no values for must fail loudly, naming the residue.
  expect_error(.aa.property.matrix("MSWHIM", c("A", "X", "Z")), "X, Z")

  # `pK` lives in Peptides::AAdata but covers only 9 residues; it used to fall
  # through as a silently wrong 9-column matrix.
  expect_error(.aa.property.matrix("pK"), "no values for")
})

test_that("built-in scales are stored in canonical amino-acid order", {
  # Guards the single-definition/collation arrangement of `amino.acids`:
  # `.builtin_scales` builds its dimnames from it at load time.
  expect_identical(colnames(.builtin_scales$atchleyFactors), amino.acids)
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
