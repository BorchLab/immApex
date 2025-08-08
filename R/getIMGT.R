#' Get IMGT Sequences for Specific Loci
#'
#' Use this to access the ImMunoGeneTics (IMGT) sequences for a
#' specific species and gene loci. More information on 
#' IMGT can be found at \href{https://www.imgt.org/}{imgt.org}.
#'
#' @examples
#' \dontrun{
#' TRBV_aa <- getIMGT(species = "human",
#'                    chain = "TRB",
#'                    frame = "inframe",
#'                    region = "v",
#'                    sequence.type = "aa", 
#'                    max.retries = 3) 
#' }
#'
#' @param species One or two-word common designation of species. 
#' @param chain Sequence chain to access, e.g., \strong{TRB} or \strong{IGH}.
#' @param frame Designation for \strong{all}, \strong{inframe}, or \strong{inframe+gap}.
#' @param region Gene loci to access.
#' @param sequence.type Type of sequence - \strong{aa} (amino acid) or \strong{nt} (nucleotide).
#' @param max.retries Number of attempts to fetch data in case of failure.
#' @param verbose Print messages corresponding to the processing step.
#'
#' @importFrom httr GET content user_agent
#' @importFrom rvest read_html html_text html_nodes
#' @importFrom stringr str_replace_all
#' @export
#' @return A list of allele sequences.
getIMGT <- function(species = "human",
                    chain = "TRB",
                    sequence.type = "aa",
                    frame = "inframe",
                    region = "v", 
                    max.retries = 3,
                    verbose = TRUE) {
  
  if (is.null(getOption("getIMGT_first_run"))) {
    message("Data from IMGT is under a CC BY-NC-ND 4.0 license. Attribution is required.")
    options(getIMGT_first_run = TRUE)
  }
  
  validate_input(region, c("v", "d", "j", "c"), "region")
  validate_input(sequence.type, c("aa", "nt"), "sequence.type")
  validate_input(frame, c("all", "inframe", "inframe+gap"), "frame")
  
  if (frame == "inframe+gap" & chain %in% c("IGLJ", "IGKJ", "IGHJ", "IGHD")) {
    stop("IMGT-gapped sequences are not available for 'IGLJ', 'IGKJ', 'IGHJ', or 'IGHD'")
  }
  
  selection <- switch(
    paste0(tolower(frame), "_", tolower(sequence.type)), 
    "all_nt" = 7.2, "inframe_nt" = 7.5, 
    "inframe_aa" = 7.6, "inframe+gap_nt" = 7.1, 
    "inframe+gap_aa" = 7.3, stop("Invalid frame and sequence type combination.")
  )
  
  chain.update <- toupper(paste0("+", chain, region))
  species.update <- str_replace_all(.parseSpecies(species), " ", "+")
  base.url <- "https://www.imgt.org/genedb/GENElect?query="
  updated.url <- paste0(base.url, selection, chain.update, "&species=", species.update)
  
  if (verbose) message("Getting sequences from IMGT...")
  
  # Attempt to fetch the webpage with retries
  response <- NULL
  attempt <- 1
  success <- FALSE
  
  while (attempt <= max.retries && !success) {
    response <- httr::GET(updated.url)
    
    if (!httr::http_error(response)) {
      success <- TRUE
    } else {
      if (verbose) {
        message(sprintf("Attempt %d failed. Retrying...", attempt))
      }
      attempt <- attempt + 1
      Sys.sleep(2)  # Optional: Add a short delay between attempts
    }
  }
  
  if (!success) {
    warning("Failed to retrieve data after ", max.retries, 
            " attempts. The website may be down or unavailable.")
    return(NULL)  
  }
  
  webpage <- httr::content(response, as = "text")
  html <- read_html(webpage)
  pre_text <- html_text(html_nodes(html, "pre"))[2]
  sequences <- parse_sequences(pre_text, chain, region, sequence.type)
  
  if (verbose) message("Formatting IMGT sequences...")
  
  list(
    sequences = sequences,
    misc = list(
      species = species.update, chain = chain,
      sequence.type = sequence.type, frame = frame, region = region
    )
  )
}

# Helper function to validate input parameters
validate_input <- function(value, valid_options, parameter_name) {
  if (tolower(value) %!in% valid_options) {
    stop(sprintf("Invalid %s. Choose one of: %s", 
                 parameter_name, paste(valid_options, collapse = ", ")))
  }
}

# Helper function to parse sequences
#' @importFrom stringr str_split str_remove_all str_extract
parse_sequences <- function(pre_text, chain, region, sequence.type) {
  sequences <- str_split(str_remove_all(pre_text, "\n"), ">")[[1]]
  fasta_list <- list()
  
  for (seq in sequences[nchar(sequences) > 0]) {
    parts <- str_split(seq, "\\|")[[1]]
    if (length(parts) < 2) next
    
    allele_name <- parts[2]
    if (tolower(sequence.type) == "aa") {
      sequence <- gsub("[^A-Z]", "", str_extract(seq, "(?<=\\|)[\\s\\S]*$"))
    } else {
      sequence <- gsub("[^acgt]", "", str_extract(seq, "[acgt]+$"))
    }
    
    fasta_list[[allele_name]] <- sequence
  }
  
  return(fasta_list)
}

#' @importFrom hash hash
.parseSpecies <- function(x) {
  species <- c("human", "mouse", "rat", "rabbit", "rhesus monkey", 
               "sheep", "pig", "platypus", "alpaca", "dog", 
               "chicken", "ferret")
  
  species_dictionary <- hash::hash(
    "human" = "Homo sapiens", 
    "mouse" = "Mus", 
    "rat" = "Rattus norvegicus", 
    "rabbit" = "Oryctolagus cuniculus", 
    "rhesus monkey" = "Macaca mulatta", 
    "sheep" = "Ovis aries", 
    "pig" = "Sus scrofa", 
    "platypus" = "Ornithorhynchus anatinus", 
    "alpaca" = "Vicugna pacos", 
    "dog" = "Canis lupus familiaris", 
    "chicken" = "Gallus gallus", 
    "ferret" = "Mustela putorius furo"
  )
  
  x <- tolower(x)
  if (x %in% species) {
    return(species_dictionary[[x]])
  } else {
    stop(sprintf("Invalid species. Choose one of: %s", 
                 paste(species, collapse = ", ")))
  }
}
