# test script for buildNetwork.R - testcases are NOT comprehensive!

# Helper function to create test data
create_test_data <- function(n = 100, 
                             n_v_genes = 5, 
                             n_j_genes = 3) {
  set.seed(42)
  
  # Generate random amino acid sequences
  aa <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
          "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  
  sequences <- replicate(n, {
    len <- sample(12:18, 1)
    paste(sample(aa, len, replace = TRUE), collapse = "")
  })
  
  # Assign V and J genes
  v_genes <- paste0("IGHV", sample(1:n_v_genes, n, replace = TRUE))
  j_genes <- paste0("IGHJ", sample(1:n_j_genes, n, replace = TRUE))
  
  data.frame(
    sequence_id = paste0("seq_", 1:n),
    junction_aa = sequences,
    v_call = v_genes,
    j_call = j_genes,
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# Test 1: Backward Compatibility
# ============================================================================

test_that("Backward compatibility - default parameters work", {
  data <- create_test_data(50)
  
  # Original usage should still work
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 2,
    filter.v = TRUE
  )
  
  expect_s3_class(edges, "data.frame")
  expect_true(all(c("from", "to", "dist") %in% names(edges)))
  expect_true(all(edges$dist <= 2))
})

test_that("Backward compatibility - relative threshold", {
  data <- create_test_data(50)
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 0.9,
    filter.v = TRUE
  )
  
  expect_s3_class(edges, "data.frame")
  expect_true(all(edges$dist >= 0 & edges$dist <= 1))
})

# ============================================================================
# Test 2: Multiple Distance Metrics
# ============================================================================

test_that("Levenshtein metric works", {
  data <- create_test_data(30)
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 2,
    dist_type = "levenshtein",
    filter.v = TRUE
  )
  
  expect_s3_class(edges, "data.frame")
  expect_true(all(edges$dist <= 2))
})

test_that("Hamming metric works with equal-length sequences", {
  # Create equal-length sequences
  data <- data.frame(
    sequence_id = paste0("seq_", 1:20),
    junction_aa = c(
      rep("ACDEFGHIKLM", 10),
      rep("ACDEFGHIKLM", 10)  # Ensure we have matches
    ),
    v_call = rep("IGHV1", 20),
    stringsAsFactors = FALSE
  )
  # Add some variations
  data$junction_aa[2] <- "ACDEFGHIKLN"  # 1 substitution
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 2,
    dist_type = "hamming",
    filter.v = TRUE
  )
  
  expect_s3_class(edges, "data.frame")
})

test_that("Damerau metric works", {
  data <- create_test_data(30)
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 2,
    dist_type = "damerau",
    filter.v = TRUE
  )
  
  expect_s3_class(edges, "data.frame")
  expect_true(all(edges$dist <= 2))
})

test_that("Needleman-Wunsch metric works", {
  data <- create_test_data(20)
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 5,
    dist_type = "nw",
    dist_mat = "BLOSUM62",
    filter.v = TRUE
  )
  
  expect_s3_class(edges, "data.frame")
})

test_that("Smith-Waterman metric works", {
  data <- create_test_data(20)
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 5,
    dist_type = "sw",
    dist_mat = "BLOSUM62",
    filter.v = TRUE
  )
  
  expect_s3_class(edges, "data.frame")
})

# ============================================================================
# Test 3: Normalization Options
# ============================================================================

test_that("No normalization (default) works", {
  data <- create_test_data(30)
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 3,
    normalize = "none",
    filter.v = TRUE
  )
  
  expect_s3_class(edges, "data.frame")
  expect_true(all(edges$dist >= 0))
  expect_true(all(edges$dist == floor(edges$dist)))  # Should be integers
})

test_that("maxlen normalization works", {
  data <- create_test_data(30)
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 0.2,
    normalize = "maxlen",
    filter.v = TRUE
  )
  
  expect_s3_class(edges, "data.frame")
  expect_true(all(edges$dist >= 0 & edges$dist <= 1))
})

test_that("length normalization works", {
  data <- create_test_data(30)
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 0.2,
    normalize = "length",
    filter.v = TRUE
  )
  
  expect_s3_class(edges, "data.frame")
  expect_true(all(edges$dist >= 0 & edges$dist <= 1))
})

# ============================================================================
# Test 4: Substitution Matrices
# ============================================================================

test_that("BLOSUM62 matrix works", {
  data <- create_test_data(20)
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 5,
    dist_type = "nw",
    dist_mat = "BLOSUM62",
    filter.v = TRUE
  )
  
  expect_s3_class(edges, "data.frame")
})

test_that("PAM30 matrix works", {
  data <- create_test_data(20)
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 5,
    dist_type = "nw",
    dist_mat  = "PAM30",
    filter.v = TRUE
  )
  
  expect_s3_class(edges, "data.frame")
})

# ============================================================================
# Test 5: V/J Filtering and Memory Efficiency
# ============================================================================

test_that("V filtering works correctly", {
  data <- create_test_data(50, n_v_genes = 3)
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 2,
    filter.v = TRUE
  )
  
  # Verify that edges only connect sequences with same V gene
  for (i in 1:nrow(edges)) {
    from_v <- data$v_call[data$sequence_id == edges$from[i]]
    to_v <- data$v_call[data$sequence_id == edges$to[i]]
    expect_equal(from_v, to_v)
  }
})

test_that("J filtering works correctly", {
  data <- create_test_data(50, n_j_genes = 3)
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 2,
    filter.j = TRUE
  )
  
  # Verify that edges only connect sequences with same J gene
  for (i in 1:nrow(edges)) {
    from_j <- data$j_call[data$sequence_id == edges$from[i]]
    to_j <- data$j_call[data$sequence_id == edges$to[i]]
    expect_equal(from_j, to_j)
  }
})

test_that("V+J filtering works correctly", {
  data <- create_test_data(50, n_v_genes = 3, n_j_genes = 2)
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 2,
    filter.v = TRUE,
    filter.j = TRUE
  )
  
  # Verify that edges only connect sequences with same V and J genes
  for (i in 1:nrow(edges)) {
    from_v <- data$v_call[data$sequence_id == edges$from[i]]
    to_v <- data$v_call[data$sequence_id == edges$to[i]]
    from_j <- data$j_call[data$sequence_id == edges$from[i]]
    to_j <- data$j_call[data$sequence_id == edges$to[i]]
    
    expect_equal(from_v, to_v)
    expect_equal(from_j, to_j)
  }
})

# ============================================================================
# Test 6: Output Formats
# ============================================================================

test_that("Edge list output works", {
  data <- create_test_data(30)
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 2,
    output = "edges",
    filter.v = TRUE
  )
  
  expect_s3_class(edges, "data.frame")
  expect_true(all(c("from", "to", "dist") %in% names(edges)))
})

test_that("Sparse matrix output works", {
  data <- create_test_data(30)
  
  mat <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 2,
    output = "sparse",
    weight = "dist",
    filter.v = TRUE
  )
  
  expect_s4_class(mat, "dgCMatrix")
})

test_that("Binary weight works with sparse matrix", {
  data <- create_test_data(30)
  
  mat <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 2,
    output = "sparse",
    weight = "binary",
    filter.v = TRUE
  )
  
  expect_s4_class(mat, "dgCMatrix")
  # All non-zero values should be 1
  expect_true(all(mat@x == 1))
})

# ============================================================================
# Test 7: Edge Cases
# ============================================================================

test_that("Empty sequences are handled", {
  data <- create_test_data(30)
  data$junction_aa[1] <- ""  # Add empty sequence
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 2,
    filter.v = TRUE
  )
  
  expect_s3_class(edges, "data.frame")
  # Empty sequence should not appear in results
  expect_false("seq_1" %in% edges$from)
  expect_false("seq_1" %in% edges$to)
})

test_that("Single sequence per group handled", {
  data <- data.frame(
    sequence_id = paste0("seq_", 1:5),
    junction_aa = c("AAAA", "BBBB", "CCCC", "DDDD", "EEEE"),
    v_call = paste0("IGHV", 1:5),  # Each sequence in different V group
    stringsAsFactors = FALSE
  )
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 2,
    filter.v = TRUE
  )
  
  # Should return empty edge list (no sequences share V gene)
  expect_s3_class(edges, "data.frame")
  expect_equal(nrow(edges), 0)
})

test_that("Very similar sequences detected", {
  data <- data.frame(
    sequence_id = paste0("seq_", 1:3),
    junction_aa = c("ACDEFGHIKLM", "ACDEFGHIKLM", "ACDEFGHIKLN"),
    v_call = rep("IGHV1", 3),
    stringsAsFactors = FALSE
  )
  
  edges <- buildNetwork(
    input.data = data,
    seq_col = "junction_aa",
    threshold = 1,
    filter.v = TRUE
  )
  
  expect_s3_class(edges, "data.frame")
  expect_true(nrow(edges) >= 2)  # Should find seq1==seq2 and seq2!=seq3
})

# ============================================================================
# Test 8: Performance and Scalability
# ============================================================================

test_that("Large dataset with V/J filtering completes", {
  skip_on_cran()  # Skip on CRAN due to time
  
  data <- create_test_data(n = 1000, n_v_genes = 20, n_j_genes = 5)
  
  # Should complete without error and use reasonable memory
  expect_silent({
    edges <- buildNetwork(
      input.data = data,
      seq_col = "junction_aa",
      threshold = 0.15,
      normalize = "maxlen",
      filter.v = TRUE,
      filter.j = TRUE
    )
  })
  
  expect_s3_class(edges, "data.frame")
})

# ============================================================================
# Test 9: Parameter Validation
# ============================================================================

test_that("Invalid normalization rejected", {
  data <- create_test_data(30)
  
  expect_error(
    buildNetwork(
      input.data = data,
      seq_col = "junction_aa",
      threshold = 2,
      normalize = "invalid_norm"
    ),
    "should be one of"
  )
})

test_that("Invalid substitution matrix rejected", {
  data <- create_test_data(20)
  
  expect_error(
    buildNetwork(
      input.data = data,
      seq_col = "junction_aa",
      threshold = 5,
      dist_type = "nw",
      dist_mat = "invalid_matrix"
    ),
    "Cannot find matrix for method"
  )
})

test_that("Hamming requires equal length", {
  data <- data.frame(
    sequence_id = c("seq1", "seq2", "seq3"),
    junction_aa = c("AAAA", "BBBB", "CCCCC"),  # Different lengths
    v_call = rep("IGHV1", 3),
    stringsAsFactors = FALSE
  )
  
  expect_warning(
    edges <- buildNetwork(
      input.data = data,
      seq_col = "junction_aa",
      threshold = 2,
      dist_type = "hamming",
      filter.v = TRUE
    ),
    "equal-length"
  )
})
