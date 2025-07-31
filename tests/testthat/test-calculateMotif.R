# test script for calculateMotif.R - testcases are NOT comprehensive!

# ---------------------------------------------------------------------------
# Basic contiguous motif counting
# ---------------------------------------------------------------------------
test_that("contiguous motifs counted correctly", {
  
  seqs <- c("AAA", "AAB", "ABA")
  res  <- calculateMotif(seqs, motif.lengths = 2, min.depth = 1,
                         discontinuous = FALSE, nthreads = 1)
  
  expect_s3_class(res, "data.frame")
  expect_named(res, c("motif", "frequency"))
  
  # For manual reference:
  # "AA": 3, "AB": 2, "BA": 1
  ref <- data.frame(
    motif     = c("AA", "AB", "BA"),
    frequency = c(3L,   2L,   1L),
    stringsAsFactors = FALSE
  )
  
  res <- res[order(res$motif), ]
  ref <- ref[order(ref$motif), ]
  
  expect_equal(res$frequency, ref$frequency)
  expect_equal(res$motif, ref$motif)
})

# ---------------------------------------------------------------------------
#  min.depth filtering works
# ---------------------------------------------------------------------------
test_that("min.depth removes low-frequency motifs", {
  
  seqs <- c("ABCABC")                  
  res  <- calculateMotif(seqs, motif.lengths = 2,
                         min.depth = 3, discontinuous = FALSE)
  
  # Nothing reaches a depth of 3
  expect_equal(nrow(res), 0)
})

# ---------------------------------------------------------------------------
#  Discontinuous (single-gap) motifs
# ---------------------------------------------------------------------------
test_that("discontinuous motifs included when requested", {
  
  seqs <- rep("CASS", 5)
  res1 <- calculateMotif(seqs, motif.lengths = 3,
                         discontinuous = FALSE)
  
  res2 <- calculateMotif(seqs, motif.lengths = 3,
                         discontinuous = TRUE, discontinuous.symbol = ".")
  
  # `res2` should contain everything in `res1` *plus* gap variants
  expect_true(all(res1$motif %in% res2$motif))
  expect_gt(nrow(res2), nrow(res1))
  
  gapped <- res2[res2$motif == "C.S", "frequency", drop = TRUE]
  expect_equal(gapped, 5L)
})

# ---------------------------------------------------------------------------
# Threading gives identical results
# ---------------------------------------------------------------------------
test_that("multi-thread and single-thread results are identical", {
  
  skip_if_not_installed("Rcpp")   # threading check irrelevant otherwise
  
  seqs <- rep(c("ARNT", "ARND"), 100)  # moderate input
  single <- calculateMotif(seqs, motif.lengths = 2:3, nthreads = 1)
  multi  <- calculateMotif(seqs, motif.lengths = 2:3, nthreads = 2)
  
  # Sort for reproducible comparison
  single <- single[order(single$motif), ]
  multi  <- multi[order(multi$motif), ]
  
  expect_equal(single, multi)
})

# ---------------------------------------------------------------------------
# Edge-cases & argument validation
# ---------------------------------------------------------------------------
test_that("argument validation works", {
  
  expect_error(calculateMotif(character(0)),           "empty")
  expect_error(calculateMotif("AAA", motif.lengths = 0),  "positive")
  expect_error(calculateMotif("AAA", min.depth = 0),      ">= 1")
  expect_error(calculateMotif("AAA", discontinuous.symbol = "XX"),
               "single character")
})

## ---------------------------------------------------------------------------
## 6.  Motif lengths longer than sequences are skipped gracefully
## ---------------------------------------------------------------------------
test_that("motif lengths exceeding sequence length drop silently", {
  
  res <- calculateMotif("ABC", motif.lengths = c(2, 4))
  # Only length-2 motifs should be present
  expect_true(all(nchar(res$motif) == 2))
})
