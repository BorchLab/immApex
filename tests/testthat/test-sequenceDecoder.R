# test script for sequenceDecoder.R - testcases are NOT comprehensive!

test_that("sequenceDecoder rejects invalid inputs", {
  # Error on invalid object type
  expect_error(
    sequenceDecoder(c(1, 2, 3)),
    "'encoded.object' must be a list from sequenceEncoder, a 3D array, or a 2D matrix."
  )
  # Error on invalid call.threshold
  expect_error(
    sequenceDecoder(array(0, dim=c(1,1,1)), call.threshold = 0),
    "'call.threshold' must be a positive number."
  )
})

test_that("sequenceDecoder handles property mode argument errors", {
  encoded <- sequenceEncoder("A", mode = "property", property.set = "atchleyFactors")
  
  # Error when no property info is provided
  expect_error(
    sequenceDecoder(encoded$cube, mode = "property"),
    "In 'property' mode, you must supply either 'property.set' or 'property.matrix'."
  )
  # Error on flattened matrix without property info
  expect_error(
    sequenceDecoder(encoded$flattened, mode = "property"),
    "For flattened matrix input in 'property' mode, supply 'property.set' or 'property.matrix'."
  )
})

test_that(".onehotDecoder decodes standard sequences correctly", {
  sequences <- c("CARD", "WYN")
  encoded <- sequenceEncoder(sequences, mode = "onehot")
  
  # Decode from the full list object
  expect_equal(sequenceDecoder(encoded, mode = "onehot"), c("CARD", "WYN"))
  
  # Decode from the raw 3D cube
  expect_equal(sequenceDecoder(encoded$cube, mode = "onehot"), c("CARD", "WYN"))
  
  # Decode from the flattened 2D matrix
  expect_equal(sequenceDecoder(encoded$flattened, mode = "onehot"), c("CARD", "WYN"))
})

test_that(".onehotDecoder handles padding and thresholds", {
  sequences <- c("CA", "R")
  encoded <- sequenceEncoder(sequences, mode = "onehot")
  
  # Test padding removal (default)
  expect_equal(sequenceDecoder(encoded$cube, mode = "onehot"), c("CA", "R"))
  
  # Test keeping padding
  expect_equal(sequenceDecoder(encoded$cube, mode = "onehot", remove.padding = FALSE), c("CA", "R."))
  
  # --- Test call.threshold ---
  ambiguous_cube <- encoded$cube
  # Create an ambiguous position 
  ambiguous_cube[c(2, 5), 1, 1] <- 0.3 
  
  # With default threshold (0.5), this is a tie and should be padded
  expect_equal(sequenceDecoder(ambiguous_cube, mode = "onehot", call.threshold = 0.5)[1], ".A")
})


test_that(".propertyDecoder decodes standard sequences correctly", {
  sequences <- c("CARD", "WYN")
  encoded <- sequenceEncoder(sequences, 
                             mode = "property", 
                             property.set = "kideraFactors")
  
  # Decode from the raw 3D cube
  expect_equal(
    as.vector(sequenceDecoder(encoded$cube, 
                    mode = "property", 
                    property.set = "kideraFactors", 
                    call.threshold = 0.1)),
    c("CARD", "WYN")
  )
  
  # Decode from the flattened 2D matrix
  expect_equal(
    sequenceDecoder(encoded$flattened, 
                    mode = "property", 
                    property.set = "kideraFactors", 
                    call.threshold = 0.1),
    c("CARD", "WYN")
  )
})

test_that(".propertyDecoder handles padding and thresholds", {
  sequences <- c("CA", "R")
  encoded <- sequenceEncoder(sequences, 
                             mode = "property", 
                             property.set = "atchleyFactors")
  
  # Test padding removal (default)
  expect_equal(
    as.vector(sequenceDecoder(encoded$cube, 
                    mode = "property", 
                    property.set = "atchleyFactors", 
                    call.threshold = 0.1)),
    c("CA", "R")
  )
  
  # Test keeping padding
  expect_equal(
    as.vector(sequenceDecoder(encoded$cube, 
                    mode = "property", 
                    property.set = "atchleyFactors", 
                    remove.padding = FALSE, 
                    call.threshold = 0.1)),
    c("CA", "R.")
  )
  
  #  Test call.threshold 
  perturbed_cube <- encoded$cube
  perturbed_cube[, 1, 1] <- perturbed_cube[, 1, 1] + 0.01 # Perturb 'C' in "CA"
  
  # With a lenient threshold, it should still be called correctly
  expect_equal(
    as.vector(sequenceDecoder(perturbed_cube, 
                    mode = "property", 
                    property.set = "atchleyFactors", 
                    call.threshold = 0.1)[1]),
    "CA"
  )
  
  # With a very strict threshold, the perturbed vector should fail and be padded
  expect_equal(
    as.vector(sequenceDecoder(perturbed_cube, 
                    mode = "property", 
                    property.set = "atchleyFactors", 
                    call.threshold = 0.001)[1]),
    ".A"
  )
})
