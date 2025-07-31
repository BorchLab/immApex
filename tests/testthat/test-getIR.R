# test script for getIR.R - testcases are NOT comprehensive!

  test_that("Test for valid sequence type", {
    dummy_data <- data.frame(barcode = c("cell1", "cell2"),
                             CTaa = c("CASSLGG", "CASRLGG"),
                             CTnt = c("TGTGCCAGCAGCTTGGG", "TGTGCCAGCAGGCTGGG"),
                             CTgene = c("TRAV1-2.TRAC", "TRBV2.TRBJ1-1"))
    
    expect_error(getIR(input.data = dummy_data, 
                       chains = "TRA",
                       sequence.type = "invalidType"),
                 regexp = "should be one of")
  })
  
  test_that("Test for valid chain type", {
    dummy_data <- data.frame(barcode = c("cell1", "cell2"),
                             CTaa = c("CASSLGG_XXXXX", "CASRLGG_XXXXX"),
                             CTnt = c("TGTGCCAGCAGCTTGGG_YYYYY", "TGTGCCAGCAGGCTGGG_YYYYY"),
                             CTgene = c("TRAV1-2.TRAJ.TRAC_ZZZZ.ZZZZZ.ZZZZ", "TRAV1-2.TRAJ.TRAC_ZZZZZ.ZZZZZ.ZZZZ"))
    
    expect_error(getIR(input.data = dummy_data, 
                       chains = "invalidChain"),
                 "`chains` must be one of: TRA, TRB, TRG, TRD, Heavy, Light")
  })
  
  test_that("Test for correct output type", {
    dummy_data <- data.frame(barcode = c("cell1", "cell2"),
                             CTaa = c("CASSLGG_XXXXX", "CASRLGG_XXXXX"),
                             CTnt = c("TGTGCCAGCAGCTTGGG_YYYYY", "TGTGCCAGCAGGCTGGG_YYYYY"),
                             CTgene = c("TRAV1-2.TRAJ.TRAC_ZZZZ.ZZZZZ.ZZZZ", "TRAV1-2.TRAJ.TRAC_ZZZZZ.ZZZZZ.ZZZZ"))
    
    result <- getIR(input.data = dummy_data, 
                    chains = "TRA")
    expect_true(is.data.frame(result))
  })
  
  test_that("Test for correct column names", {
    dummy_data <- data.frame(barcode = c("cell1", "cell2"),
                             CTaa = c("CASSLGG_XXXXX", "CASRLGG_XXXXX"),
                             CTnt = c("TGTGCCAGCAGCTTGGG_YYYYY", "TGTGCCAGCAGGCTGGG_YYYYY"),
                             CTgene = c("TRAV1-2.TRAJ.TRAC_ZZZZ.ZZZZZ.ZZZZ", "TRAV1-2.TRAJ.TRAC_ZZZZZ.ZZZZZ.ZZZZ"))
    
    result <- getIR(input.data = dummy_data, 
                    chains = "TRA")
    expected_colnames <- c("cdr3_aa", "v", "d", "j", "c", "barcode", "chain")
    expect_equal(colnames(result), expected_colnames)
  })
  
  test_that("Test for handling SingleCellExperiment object", {
    dummy_data <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = matrix(1:4, ncol = 2)),
      colData = data.frame(barcode = c("cell1", "cell2"),
                           CTaa = c("CASSLGG_XXXXX", "CASRLGG_XXXXX"),
                           CTnt = c("TGTGCCAGCAGCTTGGG_YYYYY", "TGTGCCAGCAGGCTGGG_YYYYY"),
                           CTgene = c("TRAV1-2.TRAJ.TRAC_ZZZZ.ZZZZZ.ZZZZ", "TRAV1-2.TRAJ.TRAC_ZZZZZ.ZZZZZ.ZZZZ"))
    )
    
    result <- getIR(input.data = dummy_data, 
                    chains = "TRA")
    expect_true(is.data.frame(result))
  })
  
