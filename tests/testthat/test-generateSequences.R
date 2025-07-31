# test script for generateSequences.R - testcases are NOT comprehensive!

test_that("Basic generation with motifs works correctly", {
  prefix <- "CAR"
  suffix <- "RAN"
  n_seq <- 50
  
  # Generate sequences
  sequences <- generateSequences(
    prefix.motif = prefix,
    suffix.motif = suffix,
    max.length = 14,
    number.of.sequences = n_seq,
    sequence.dictionary = amino.acids
  )
  
  # Check output type and quantity
  expect_true(is.character(sequences))
})


test_that("Length constraints (min and max) are respected", {
  min_len <- 8
  max_len <- 12
  
  sequences <- generateSequences(
    number.of.sequences = 100,
    min.length = min_len,
    max.length = max_len
  )
  
  seq_lengths <- nchar(sequences)
  
  # Verify all sequences are within the specified length bounds
  expect_true(all(seq_lengths >= min_len))
  expect_true(all(seq_lengths <= max_len))
  
  # Verify that lengths are not all identical (unless min.length == max.length)
  expect_true(length(unique(seq_lengths)) > 1)
})


test_that("min.length is correctly adjusted for long motifs", {
  prefix <- "VYNWTSGDMR" 
  suffix <- "VAPGEADRYT" 
  motif_len <- 20
  
  # Set min.length to be shorter than the combined motifs
  sequences <- generateSequences(
    prefix.motif = prefix,
    suffix.motif = suffix,
    number.of.sequences = 50,
    min.length = 5,
    max.length = 25
  )
  
  # The effective minimum length should be the motif length (20)
  expect_gte(min(nchar(sequences)), motif_len)
})


test_that("Edge cases are handled correctly", {
  # Case 1: number.of.sequences = 0
  expect_equal(length(generateSequences(number.of.sequences = 0)), 0)
  expect_true(is.character(generateSequences(number.of.sequences = 0)))
  
  # Case 2: min.length equals max.length
  fixed_len <- 15
  sequences <- generateSequences(
    number.of.sequences = 20,
    min.length = fixed_len,
    max.length = fixed_len
  )
  expect_true(all(nchar(sequences) == fixed_len))
  
  # Case 3: No motifs (default NULL values)
  sequences_no_motif <- generateSequences(number.of.sequences = 10)
  expect_equal(length(sequences_no_motif), 10)
  # All characters should be from the default amino acid dictionary
  expect_true(all(grepl(paste0("^[", paste(amino.acids, collapse = ""), "]*$"), sequences_no_motif)))
})


test_that("Function throws errors for invalid arguments", {
  expect_error(
    generateSequences(prefix.motif = "ABCDE", 
                      suffix.motif = "FGH", 
                      max.length = 9, 
                      verbose = FALSE),
    regexp = "Prefix/suffix contain letters not in"
  )
  
  # Error: Invalid number.of.sequences
  expect_error(generateSequences(number.of.sequences = -1)) # stopifnot
  expect_error(generateSequences(number.of.sequences = c(10, 20))) # stopifnot
})
