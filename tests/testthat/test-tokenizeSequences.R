# test script for tokenizeSequences.R - testcases are NOT comprehensive!

aa <- amino.acids
seqs <- c("ACD",            # 3-mer
          "FGHIK")          # 5-mer

test_that("matrix mode returns correct type, dimensions, and padding", {
  mat <- tokenizeSequences(seqs, 
                           add.startstop = TRUE, 
                           convert.to.matrix = TRUE,
                           verbose = FALSE)
  
  expect_true(is.matrix(mat) && is.integer(mat))
  
  # rows == sequences, cols == longest(seq)+2 (start+stop)
  expect_equal(dim(mat), c(length(seqs), max(nchar(seqs)) + 2))
  
  # default pad id = |char_set| + 1
  pad.id <- length(c("!", aa, "^")) + 1L
  # last two cells of first row must be padding
  expect_equal(mat[1, (nchar(seqs[1]) + 3):ncol(mat)], rep(pad.id, 2))
})

test_that("list mode mirrors matrix logic", {
  lst <- tokenizeSequences(seqs,
                           add.startstop = FALSE,
                           convert.to.matrix = FALSE,
                           max.length = 6,
                           padding.symbol = 99,
                           verbose = FALSE)
  expect_type(lst, "list")
  expect_length(lst, length(seqs))
  expect_equal(lengths(lst), rep(6L, 2))          # padded to 6
})

test_that("max.length shorter than longest sequence triggers error", {
  expect_error(
    tokenizeSequences(seqs, 
                      max.length = 2,
                      verbose = FALSE),
    "`max.length` is smaller"
  )
})

test_that("unknown characters are caught", {
  bad <- c("AXZ")  # Z not in amino.acids
  expect_error(
    tokenizeSequences(bad, 
                      add.startstop = FALSE,
                      verbose = FALSE),
    "Unknown character"
  )
})

test_that("non-single-character tokens are rejected", {
  expect_error(
    tokenizeSequences("AC", 
                      start.token = "XX",
                      verbose = FALSE),
    "single characters"
  )
  expect_error(
    tokenizeSequences("AC", 
                      stop.token  = "YY",
                      verbose = FALSE),
    "single characters"
  )
})

test_that("empty input returns empty structure of appropriate type", {
  m <- tokenizeSequences(character(0))
  expect_true(is.matrix(m) && length(m) == 0)
  
  l <- tokenizeSequences(character(0), convert.to.matrix = FALSE)
  expect_true(is.list(l) && length(l) == 0)
})

test_that("custom padding symbol honoured in matrix output", {
  mat <- tokenizeSequences(seqs,
                           convert.to.matrix = TRUE,
                           padding.symbol = 777)
  expect_true(all(mat[mat == 777L] == 777L))
})

test_that("all returned integers map to the declared character set", {
  mat <- tokenizeSequences(seqs, convert.to.matrix = TRUE)
  pad.id <- max(mat)
  # ids excluding padding
  ids <- setdiff(unique(as.vector(mat)), pad.id)
  expect_true(all(ids %in% seq_along(c("!", aa, "^"))))
})

