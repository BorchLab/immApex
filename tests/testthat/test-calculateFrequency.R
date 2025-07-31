# test script for calculateFrequency.R - testcases are NOT comprehensive!

## Canonical residue sets ---------------------------------------------------
aa_20 <- c("A","C","D","E","F","G","H","I","K","L",
           "M","N","P","Q","R","S","T","V","W","Y")
dna_4 <- c("A","C","G","T")

# Helper to check column sums ≈ 1 -----------------------------------------
cols_sum_to_one <- function(mat, tol = 1e-12) {
  all(abs(colSums(mat) - 1) < tol)
}

## -------------------------------------------------------------------------
test_that("matrix dimensions & names are correct (AA default)", {
  seqs <- c("CASSLGQGAETQYF", "CASSPGQGDYEQYF", "CASSQETQYF")
  res  <- calculateFrequency(seqs)
  
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(length(aa_20)+1, max(nchar(seqs))))
  expect_contains(rownames(res), aa_20)
  expect_true(cols_sum_to_one(res))
})

## -------------------------------------------------------------------------
test_that("works with nucleotide alphabet & custom padding", {
  dna <- c("ATGCC", "ATGAC", "ATGGC")
  pad <- "-"
  res <- calculateFrequency(dna,
                             sequence.dictionary = dna_4,
                             padding.symbol = pad)
  
  expect_equal(dim(res), c(5, 5))
  expect_setequal(rownames(res), c(dna_4, pad))
  expect_identical(as.vector(colSums(res)), rep(1, 5))
})

## -------------------------------------------------------------------------
test_that("`tidy = TRUE` agrees with matrix output", {
  seqs <- c("AAA", "AAC")
  mat  <- calculateFrequency(seqs, max.length = 3)
  tidy <- calculateFrequency(seqs, max.length = 3, tidy = TRUE)
  
  # reshape matrix for comparison
  mat_long <- as.data.frame(as.table(mat),
                            stringsAsFactors = FALSE,
                            responseName = "frequency")
  names(mat_long) <- c("residue", "position", "frequency")
  mat_long$position <- as.integer(sub("Pos\\.", "", mat_long$position))
  
  expect_equal(
    tidy[order(tidy$residue, tidy$position), ],
    mat_long[order(mat_long$residue, mat_long$position), ],
    tolerance = 1e-12
  )
})


## -------------------------------------------------------------------------
test_that("non-character input triggers error", {
  expect_error(calculateFrequency(1:5), "is.character")
})

test_that("padding symbol collision detected", {
  expect_error(
    calculateFrequency(c("AA"), padding.symbol = "A"),
    "padding.symbol"
  )
})

## -------------------------------------------------------------------------
test_that("unknown residues are ignored but columns remain normalised", {
  seqs <- c("XZA", "AAA")
  res  <- calculateFrequency(seqs, max.length = 3)
  
  # Unknowns (X,Y,Z) should contribute zero rows (already absent)
  expect_false(any(c("X", "Z") %in% rownames(res)))
  
  # Because unknowns appear, column sums should be < 1
  expect_true(all(colSums(res) <=1))
})

## -------------------------------------------------------------------------
test_that("single-sequence edge case returns sensible output", {
  seqs <- "CASSQETQYF"
  res  <- calculateFrequency(seqs)
  
  expect_true(is.matrix(res))
  expect_true(cols_sum_to_one(res))
})
