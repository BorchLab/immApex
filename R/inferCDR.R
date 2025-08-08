#' Infer CDR-loop segments from V-gene calls
#' 
#' Use this isolate sequences from the CDR loop using
#' the V gene annotation. When there are multiple 
#' V gene matches for a single gene, the first allelic 
#' sequence is used.
#' 
#' @examples
#'  \dontrun{
#' # Getting the Sequence Reference
#' data(immapex_example.data)
#' TRBV_aa <- getIMGT(species = "human",
#'                    chain = "TRB",
#'                    frame = "inframe",
#'                    region = "v",
#'                    sequence.type = "aa") 
#'                     
#' # Ensuring sequences are formatted to IMGT                   
#' TenX_formatted <- formatGenes(immapex_example.data[["TenX"]],
#'                               region = "v",
#'                               technology = "TenX")
#'              
#' # Inferring CDR loop elements           
#' TenX_formatted <- inferCDR(TenX_formatted,
#'                            chain = "TRB", 
#'                            reference = TRBV_aa,
#'                            technology = "TenX", 
#'                            sequence.type = "aa",
#'                            sequences = c("CDR1", "CDR2"))
#' }
#'                             
#' @param input.data Data frame output of \code{\link{formatGenes}}
#' @param reference IMGT reference sequences from \code{\link{getIMGT}}
#' @param chain Sequence chain to access, like \strong{TRB} or \strong{IGH}
#' @param technology The sequencing technology employed - \strong{TenX}, 
#' \strong{Adaptive}, or \strong{AIRR}
#' @param sequence.type Type of sequence - \strong{aa} for amino acid or 
#' \strong{nt} for nucleotide
#' @param sequences The specific regions of the CDR loop to get from the data, 
#' such as \strong{CDR1}. 
#' @param verbose Logical. If `TRUE` (default), prints a progress message.
#'                             
#' @importFrom stringr str_split
#' 
#' @export
#' @return A data frame with the new columns of CDR sequences added.
inferCDR <- function(input.data,
                     reference,
                     chain       = "TRB",
                     technology  = c("TenX", "AIRR", "Adaptive", "Omniscope"),
                     sequence.type = c("aa", "nt"),
                     sequences   = c("CDR1", "CDR2"),
                     verbose     = TRUE) {
  
  # Preflight checks-----------------------------------------------------------
  sequence.type <- match.arg(sequence.type)
  technology    <- match.arg(technology)
  
  if (is.null(reference) ||
      !"v" %in% reference[["misc"]][["region"]])
    stop("`reference` must be a V-gene list from getIMGT().")
  
  if (reference[["misc"]][["sequence.type"]] != sequence.type)
    stop("`sequence.type` mismatch with supplied reference.")
  
  unknown_regions <- setdiff(sequences, names(.sequence.positions))
  if (length(unknown_regions))
    stop("Unknown regions: ", paste(unknown_regions, collapse = ", "))
  
  # Determine V-gene column ---------------------------------------------------
  v.col <- .get_v_column(input.data, technology)
  
  # Build region positions ----------------------------------------------------
  pos_idx <- .sequence.positions[sequences]
  if (sequence.type == "nt") {
    pos_idx <- lapply(pos_idx, function(x) {
      aa_start <- min(x)
      aa_end <- max(x)
      nt_start <- (aa_start - 1L) * 3L + 1L
      nt_end <- aa_end * 3L
      seq.int(nt_start, nt_end)
    })
  }
  
  # Build a lookup table ------------------------------------------------------
  if (verbose) message("Building V-gene / CDR map...")
  
  ref_seqs <- reference[["sequences"]]
  ref_names <- names(ref_seqs)
  
  extract_one <- function(seq, idx) substring(seq, min(idx), max(idx))
  
  cdr_map <- lapply(pos_idx, function(idx) {
    vapply(ref_seqs, extract_one, idx = idx, FUN.VALUE = character(1L))
  })
  cdr_df <- data.frame(v_IMGT = ref_names,
                       do.call(cbind, cdr_map),
                       row.names = NULL,
                       check.names = FALSE,
                       stringsAsFactors = FALSE)
  names(cdr_df)[-1] <- paste0(sequences, "_IMGT")
  
  # Splice into input.data via match ------------------------------------------
  key <- .match.gene(input.data[[v.col]], cdr_df$v_IMGT)

  for (j in seq_along(sequences)) {
    new_col <- paste0(sequences[j], "_IMGT")
    input.data[[new_col]] <- cdr_df[[new_col]][key]
  }
  
  if (verbose) {
    miss <- sum(is.na(key))
    if (miss)
      message("Warning:", miss, "V genes not found in reference; CDRs set to NA.\n")
  }
  
  input.data
}


#CDR3 AA sequences by IMGT
.sequence.positions <- list(
  FR1.A = c(1:15), 
  FR1.B = c(16:26),
  CDR1 = c(27:38), 
  FR2.C = c(39:46), 
  FR2.Cp = c(47:55),
  CDR2 = c(56:65), 
  FR3.C = c(66:74), 
  FR3.D = c(75:84),
  FR3.E = c(85:96), 
  FR3.F = c(97:104)
)
