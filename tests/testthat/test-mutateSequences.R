# test script for mutateSequences.R - testcases are NOT comprehensive!

  aa <- amino.acids          
  seq1 <- "ACDEFGHIKLMNPQRSTVWY"             
  seq2 <- "CASSLGQGAETQYF"                  
  input <- c(seq1, seq2)
  
  test_that("basic output shape is correct", {
    out <- mutateSequences(input, number.of.sequences = 3, mutation.rate = 0.1)
    expect_type(out, "character")
    expect_length(out, length(input) * 3)
  })
  
  test_that("identical length is maintained", {
    out <- mutateSequences(input, number.of.sequences = 2, mutation.rate = 0.2)
    expect_equal(as.numeric(nchar(out)), rep(nchar(rep(input, each = 2)), 1))
  })
  
  test_that("mutation.rate = 0 returns originals verbatim", {
    expect_identical(
      as.character(mutateSequences(input, mutation.rate = 0)),
      input
    )
  })
  
  test_that("mutation.rate = 1 changes *every* position", {
    set.seed(42)
    mut <- mutateSequences(seq1, mutation.rate = 1)  # single string
    expect_false(any(strsplit(mut, "")[[1]] ==
                       strsplit(seq1, "")[[1]]))
  })
  
  test_that("requested window boundaries are honoured", {
    set.seed(1)
    win <- mutateSequences(seq1,
                           mutation.rate = 1,
                           position.start = 5,
                           position.end   = 10)
    unchanged <-  setdiff(seq_len(nchar(seq1)), 5:10)
    expect_equal(as.character(substr(win, unchanged, unchanged)),
                 substr(seq1, unchanged, unchanged))
  })
  
  test_that("ceil((p1-p0+1)*rate) mutations applied (probabilistic)", {
    # deterministic check by seeding and forcing small window
    set.seed(123)
    win_len <- 10
    rate    <- 0.25
    out <- mutateSequences(seq1,
                           mutation.rate = rate,
                           position.start = 1,
                           position.end   = win_len)
    num_mut <- sum(strsplit(out, "")[[1]][1:win_len] !=
                     strsplit(seq1, "")[[1]][1:win_len])
    expect_equal(num_mut, ceiling(win_len * rate))
  })
  
  test_that("all output symbols come from sequence.dictionary", {
    out <- mutateSequences(seq1, mutation.rate = 1)
    expect_true(all(unique(unlist(strsplit(out, ""))) %in% aa))
  })
  
  test_that("invalid input raises informative errors", {
    expect_error(
      mutateSequences(123),                             # non-character
      "is.character"
    )
    expect_error(
      mutateSequences(seq1, mutation.rate = -0.1),      # negative rate
      "mutation.rate"
    )
    expect_error(
      mutateSequences(seq1, position.start = 10,
                      position.end = 5),                # reversed window
      "position.start > position.end"
    )
    expect_error(
      mutateSequences("AX", sequence.dictionary = "A"), # invalid letter
      "not in"
    )
  })
  
  test_that("number.of.sequences = 0 gives zero-length result", {
    expect_length(mutateSequences(seq1, number.of.sequences = 0), 0)
  })
  
  test_that("custom dictionary is accepted and used", {
    custom <- c("A", "B", "C")
    set.seed(7)
    out <- mutateSequences("AAA", mutation.rate = 1,
                           sequence.dictionary = custom)
    expect_true(all(strsplit(out, "")[[1]] %in% custom))
  })

