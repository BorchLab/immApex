# test script for adjacencyMatrix.R - testcases are NOT comprehensive!

test_that("Correctly calculates adjacency for simple cases", {
  sequences <- c("AR", "RA")
  
  # Test with directed = TRUE and normalize = FALSE
  adj_matrix_directed <- adjacencyMatrix(sequences, 
                                         normalize = FALSE, 
                                         sequence.dictionary = c("A", "R"), 
                                         directed = TRUE)
  expected_directed <- matrix(c(0, 1, 1, 0), nrow = 2, dimnames = list(c("A", "R"), c("A", "R")))
  expect_equal(adj_matrix_directed, expected_directed)
  
  # Test with directed = FALSE (default) and normalize = FALSE
  adj_matrix_undirected <- adjacencyMatrix(sequences, 
                                           normalize = FALSE, 
                                           sequence.dictionary = c("A", "R"))
  expected_undirected <- matrix(c(0, 2, 2, 0), nrow = 2, dimnames = list(c("A", "R"), c("A", "R")))
  expect_equal(adj_matrix_undirected, expected_undirected)
})

test_that("Normalization works correctly", {
  sequences <- c("AG", "GA")
  
  # When normalized, the sum of the matrix should be 1.0
  adj_matrix_normalized <- adjacencyMatrix(sequences, normalize = TRUE, sequence.dictionary = c("A", "G"))
  expect_equal(sum(adj_matrix_normalized), 1.0)
  
  expected_normalized <- matrix(c(0, 0.5, 0.5, 0), nrow = 2, dimnames = list(c("A", "G"), c("A", "G")))
  expect_equal(adj_matrix_normalized, expected_normalized)
})

test_that("Works with a custom dictionary", {
  sequences <- c("AT", "GC")
  dna_dict <- c("A", "T", "G", "C")
  adj_matrix <- adjacencyMatrix(sequences, normalize = FALSE, sequence.dictionary = dna_dict, directed = TRUE)
  
  # Expected matrix should have dimensions 4x4
  expect_equal(dim(adj_matrix), c(4, 4))
  # Check for correct dimnames
  expect_equal(rownames(adj_matrix), dna_dict)
  expect_equal(colnames(adj_matrix), dna_dict)
  # Check values
  expect_equal(adj_matrix["A", "T"], 1)
  expect_equal(adj_matrix["G", "C"], 1)
  expect_equal(sum(adj_matrix), 2)
})


test_that("Throws errors for invalid inputs", {
  # Error for sequences with unknown characters
  sequences_bad_char <- c("AXG")
  expect_error(adjacencyMatrix(sequences_bad_char, sequence.dictionary = c("A", "G")), 
               "Unknown letters found: X")
  
  # Error for sequences that are too short to have adjacencies
  sequences_too_short <- c("A", "G")
  expect_error(adjacencyMatrix(sequences_too_short, sequence.dictionary = c("A", "G")),
               "All sequences have length < 2 - no adjacencies to compute.")
  
  # Stopifnot checks for input types
  expect_error(adjacencyMatrix(input.sequences = 123))
  expect_error(adjacencyMatrix(input.sequences = "AG", sequence.dictionary = 123))
})

test_that("Handles empty input gracefully", {
  empty_sequences <- character(0)
  adj_matrix <- adjacencyMatrix(empty_sequences, sequence.dictionary = c("A", "B"))
  
  # Should return a zero matrix with correct dimensions
  expected_matrix <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("A", "B"), c("A", "B")))
  expect_equal(adj_matrix, expected_matrix)
})

test_that("Correctly processes a more complex list of sequences", {
  sequences <- c("CASA", "CAS", "CSA")
  dict <- c("C", "A", "S")
  
  # Directed and not normalized
  adj_matrix <- adjacencyMatrix(sequences, 
                                normalize = FALSE, 
                                sequence.dictionary = dict, 
                                directed = TRUE)
  expected_matrix <- matrix(c(0, 2, 1, 0, 0, 2, 0, 2, 0), nrow = 3, byrow = TRUE,
                            dimnames = list(dict, dict))
  colnames(expected_matrix) <- dict # Add column names to expected matrix
  
  expect_equal(adj_matrix, expected_matrix)
})


test_that("Works with the default amino acid dictionary", {
  sequences <- c("VWY", "YWF")
  adj_matrix <- adjacencyMatrix(sequences, normalize = FALSE, sequence.dictionary = amino.acids, directed = TRUE)
  
  expect_equal(dim(adj_matrix), c(20, 20))
  expect_equal(adj_matrix["V", "W"], 1)
  expect_equal(adj_matrix["W", "Y"], 1)
  expect_equal(adj_matrix["Y", "W"], 1)
  expect_equal(adj_matrix["W", "F"], 1)
  # Sum of all non-zero elements
  expect_equal(sum(adj_matrix), 4)
})

