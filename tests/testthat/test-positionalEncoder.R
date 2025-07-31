# test script for positionalEncoder.R - testcases are NOT comprehensive!

test_that("Output has correct dimensions and type", {
  ml <- 50
  dm <- 64
  enc <- positionalEncoder(max.length = ml, d.model = dm)
  
  # Check type
  expect_true(is.matrix(enc))
  expect_true(is.numeric(enc))
  
  # Check dimensions
  expect_equal(dim(enc), c(ml, dm))
  expect_equal(nrow(enc), ml)
  expect_equal(ncol(enc), dm)
})

test_that("Encoding values are mathematically correct", {
  ml <- 2
  dm <- 4
  b <- 10000
  
  expected <- matrix(
    c(sin(1), cos(1), sin(0.01), cos(0.01),
      sin(2), cos(2), sin(0.02), cos(0.02)),
    nrow = 2,
    ncol = 4,
    byrow = TRUE
  )
  
  actual <- positionalEncoder(max.length = ml, d.model = dm, base = b)
  
  expect_equal(actual, expected, tolerance = 1e-9)
})


test_that("sequences argument correctly determines max.length", {
  my_seqs <- c("A", "SEQ", "LONGEST") # max length is 7
  enc <- positionalEncoder(input.sequences = my_seqs, d.model = 16)
  
  expect_equal(nrow(enc), 7)
})

test_that("explicit max.length takes precedence over sequences argument", {
  my_seqs <- c("SHORT", "SEQ") # max length is 5
  enc <- positionalEncoder(input.sequences = my_seqs, max.length = 10, d.model = 16)
  
  expect_equal(nrow(enc), 10)
})

test_that("position_offset argument works correctly", {
  enc_offset_0 <- positionalEncoder(max.length = 5, 
                                    d.model = 8, 
                                    position.offset = 0)
  
  expected_first_row <- rep(c(0, 1), times = 8 / 2)
  
  expect_equal(enc_offset_0[1, ], expected_first_row, tolerance = 1e-9)
})

test_that("base argument correctly changes the output", {
  enc_base_10k <- positionalEncoder(max.length = 10, 
                                    d.model = 16, 
                                    base = 10000)
  enc_base_100 <- positionalEncoder(max.length = 10, 
                                    d.model = 16, 
                                    base = 100)
  
  # The two matrices should be completely different
  expect_false(isTRUE(all.equal(enc_base_10k, enc_base_100)))
})


test_that("Function throws errors for invalid or missing arguments", {
  # Error: d.model is missing
  expect_error(
    positionalEncoder(max.length = 10),
    regexp = "`d.model` \\(embedding dimension\\) must be specified."
  )
  
  # Error: d.model is an odd number
  expect_error(
    positionalEncoder(max.length = 10, d.model = 7),
    regexp = "`d.model` must be an even number."
  )
  
  # Error: Neither max.length nor sequences are provided
  expect_error(
    positionalEncoder(d.model = 8),
    regexp = "Must provide either `max.length` or a set of `sequences`."
  )
})