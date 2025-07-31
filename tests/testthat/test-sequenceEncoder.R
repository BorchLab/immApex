# test script for sequenceEncoder .R - testcases are NOT comprehensive!

test_that("one-hot encoding returns correct shape and content", {
  
  seqs <- c("ACD",    # len 3
            "W",      # len 1  -> will be padded
            "MNPQR")  # len 5  -> will be truncated to max.length = 4
  
  res <- sequenceEncoder(seqs,
                         mode           = "onehot",
                         max.length     = 4,       # force trunc/pad
                         padding.symbol = ".",
                         verbose        = FALSE)
  
  expect_type(res,  "list")
  # Test for new 'mode' output and updated C++ return names
  expect_named(res, c("cube", "flattened", "sequence.dictionary", "padding.symbol", "threads", "mode"),
               ignore.order = TRUE)
  
  # depth should be alphabet + padding = 21
  expect_equal(dim(res$cube), c(21, 4, length(seqs)))
  expect_equal(dim(res$flattened), c(length(seqs), 21 * 4))
  
  ## every position has *exactly one* 1 in the depth dimension
  one.count <- as.vector(apply(res$cube, 3, function(slice) sum(slice[1:20,])))
  expect_equal(one.count, c(3,1,4))
  
  ## wrapper produces identical result
  res.w <- onehotEncoder(seqs, max.length = 4, verbose = FALSE, padding.symbol = ".")
  expect_equal(res$flattened, res.w$flattened)
})

test_that("property encoding with explicit matrix works and returns summary", {
  
  # toy property matrix: 20 AAs × 2 scales
  set.seed(42)
  amino.acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                   "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  prop.mat <- matrix(runif(40), nrow = 20,
                     dimnames = list(amino.acids, c("scale1", "scale2")))
  
  seqs <- c("AAAA", "CC", "EFG")      
  res <- sequenceEncoder(seqs,
                         mode            = "property",
                         property.matrix = prop.mat,
                         summary.fun     = "mean", 
                         verbose         = FALSE)
  
  expect_equal(dim(res$cube),       c(2, max(nchar(seqs)), length(seqs)))
  expect_equal(dim(res$flattened),  c(length(seqs), 2 * max(nchar(seqs))))
  expect_equal(dim(res$summary),    c(length(seqs), 2))
  
  res.w <- propertyEncoder(seqs,
                           property.matrix = prop.mat,
                           summary.fun     = "mean",
                           verbose         = FALSE)
  expect_equal(res$flattened, res.w$flattened)
})

test_that("geometric encoding returns correct shape and content", {
  
  seqs <- c("CARDRST", "YYYGMD", "ACACACAC")
  
  res <- sequenceEncoder(seqs,
                         mode    = "geometric",
                         method  = "BLOSUM62",
                         verbose = FALSE)
  
  # Check output structure: no cube/flattened, but has summary
  expect_type(res, "list")
  expect_named(res, c("cube", "flattened", "summary", "mode", "method"))
  expect_null(res$cube)
  expect_null(res$flattened)
  
  # Check dimensions of the summary matrix (primary output)
  expect_true(is.matrix(res$summary))
  expect_equal(dim(res$summary), c(length(seqs), 20))
  
  # Check that the alias works and gives identical output
  res.w <- geometricEncoder(seqs, method = "BLOSUM62", verbose = FALSE)
  expect_equal(res$summary, res.w$summary)
  
  # Test logical property: average of "A" should be same as "AA"
  res_A <- geometricEncoder("A", verbose = FALSE)
  res_AA <- geometricEncoder("AA", verbose = FALSE)
  expect_equal(res_A$summary, res_AA$summary)
})

test_that("geometric encoding throws informative errors", {
  
  # Error on non-canonical amino acids
  expect_error(
    sequenceEncoder("ACDX", mode = "geometric", verbose = FALSE),
    regexp = "Non-canonical"
  )
  
  # Error on invalid substitution matrix name
  expect_error(
    sequenceEncoder("ACD", mode = "geometric", method = "BOGUS_MATRIX", verbose = FALSE),
    regexp = "Cannot find matrix for method 'BOGUS_MATRIX'"
  )
})


test_that("informative errors are thrown for invalid input", {
  
  ## property mode without matrix or set
  expect_error(
    sequenceEncoder("AC", mode = "property", verbose = FALSE),
    regexp = "mode, supply either"
  )
  
  ## property.matrix wrong dimensions
  bad.mat <- matrix(1, nrow = 5, ncol = 1)
  expect_error(
    sequenceEncoder("AC",
                    mode = "property",
                    property.matrix = bad.mat,
                    verbose = FALSE),
    "rows"
  )
  
  ## zero sequences (only errors in C++ modes)
  expect_error(sequenceEncoder(character(0), mode = "onehot", verbose = FALSE),
               "`sequences` is empty", fixed = TRUE)
  
  expect_no_error(sequenceEncoder(character(0), mode = "geometric", verbose = FALSE))
})


test_that("meta slots are coherent", {
  
  seqs <- c("AC", "DFG")
  res  <- sequenceEncoder(seqs, verbose = FALSE, mode = "onehot") # specify mode
  
  # Adjusting for C++ return value names
  expect_equal(res$padding.symbol, ".")
  expect_false(res$padding.symbol %in% res$sequence.dictionary)
  expect_true(res$threads >= 1)
  expect_equal(res$mode, "onehot")
})
