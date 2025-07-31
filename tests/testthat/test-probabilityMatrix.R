# test script for probabilityMatrix .R - testcases are NOT comprehensive!


dna_sequences <- c("ATGC", "AAGC", "ATGG")
dna_dict <- c("A", "T", "G", "C")
protein_sequences <- c("CAS", "CSA")



test_that("PPM is calculated correctly for equal-length sequences", {
  ppm <- probabilityMatrix(dna_sequences, sequence.dictionary = dna_dict, convert.PWM = FALSE)
  
  # Check dimensions and names
  expect_equal(nrow(ppm), 4)
  expect_equal(ncol(ppm), 4)
  expect_equal(rownames(ppm), dna_dict)
  
  # Check that columns sum to 1
  expect_true(all(abs(colSums(ppm) - 1) < 1e-9))
  
  expect_equal(ppm["A", "Pos.1"], 1.0)
  expect_equal(ppm["T", "Pos.2"], 2/3)
  expect_equal(ppm["G", "Pos.3"], 1.0)
  expect_equal(ppm["C", "Pos.4"], 2/3)
})

test_that("PPM handles varying-length sequences and padding correctly", {
  sequences <- c("AG", "A", "AGC")
  ppm <- probabilityMatrix(sequences, sequence.dictionary = c("A","G","C"), convert.PWM = FALSE)
  
  # max.length should be 3
  expect_equal(ncol(ppm), 3)
  
  # Check column totals still sum to 1
  expect_true(all(abs(colSums(ppm) - 1) < 1e-9))
  expect_equal(ppm["A", "Pos.1"], 1.0)
  expect_equal(ppm["G", "Pos.2"], 1.0)
  expect_equal(ppm["C", "Pos.3"], 1.0)
})

test_that("PWM conversion works with uniform background", {
  sequences <- c("AA", "AC")
  dict <- c("A", "C")
  pwm <- probabilityMatrix(sequences, sequence.dictionary = dict, convert.PWM = TRUE, pseudocount = 1)
  
  expect_equal(pwm["A", "Pos.1"], log2(1.5))
  expect_equal(pwm["C", "Pos.1"], log2(0.5))
  expect_equal(pwm["A", "Pos.2"], 0)
  expect_equal(pwm["C", "Pos.2"], 0)
})

test_that("PWM conversion works with custom background frequencies", {
  sequences <- c("AA", "AC")
  dict <- c("A", "C")
  bg <- c(A = 0.8, C = 0.2)
  pwm <- probabilityMatrix(sequences, sequence.dictionary = dict, convert.PWM = TRUE, 
                           background.frequencies = bg, pseudocount = 1)
  
  expect_equal(pwm["A", "Pos.1"], log2(0.75 / 0.8))
  expect_equal(pwm["C", "Pos.1"], log2(0.25 / 0.2))
  expect_equal(pwm["A", "Pos.2"], log2(0.5 / 0.8))
  expect_equal(pwm["C", "Pos.2"], log2(0.5 / 0.2))
})

test_that("Function handles edge cases correctly", {
  # Empty input
  expect_equal(ncol(probabilityMatrix(character(0))), 0)
  
  # Column with only padding symbols
  sequences <- c("A", "C")
  ppm <- probabilityMatrix(sequences, max.length = 2, sequence.dictionary = c("A", "C"), convert.PWM = FALSE)
  expect_equal(sum(ppm[, "Pos.2"]), 0)
  
  # PWM for a column with only padding
  pwm <- probabilityMatrix(sequences, max.length = 2, sequence.dictionary = c("A", "C"), convert.PWM = TRUE)
  expect_equal(sum(abs(pwm[, "Pos.2"])), 0)
})

test_that("Function throws appropriate errors", {
  # Padding symbol is in the dictionary
  expect_error(
    probabilityMatrix("A", padding.symbol = "A", sequence.dictionary = c("A","C")),
    "`padding.symbol` cannot be present in `sequence.dictionary`."
  )
  
  # Background frequencies don't match dictionary
  expect_error(
    probabilityMatrix("AC", convert.PWM = TRUE, background.frequencies = c(A = 1)),
    "Names of `background.frequencies` must match `sequence.dictionary`."
  )
})

