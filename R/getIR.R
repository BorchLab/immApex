#' Extract Immune Receptor Sequences
#' 
#' Use this to extract immune receptor sequences from a Single-Cell 
#' Object or the output of \link[scRepertoire]{combineTCR} and 
#' \link[scRepertoire]{combineBCR}.
#' 
#' @param input.data Single-cell object or the output of 
#' \link[scRepertoire]{combineTCR} and \link[scRepertoire]{combineBCR} from
#' scRepertoire
#' @param chains Immune Receptor chain to use - \strong{TRA}, 
#' \strong{TRB}, \strong{IGH}, or \strong{IGL}
#' @param sequence.type Extract amino acid (\strong{aa}) or 
#' nucleotide (\strong{nt}) sequences
#' @param group.by Optional metadata column (e.g., \code{"sample.id"}) to 
#' group and return results as a named list by that variable.
#' @param as.list Logical; if \code{TRUE}, returns a list split by chain. 
#' If \code{group.by} is also provided, returns a nested list Default is 
#' \code{FALSE}.
#' 
#' @export
#' @return A data frame, list of data frames, or nested list of immune 
#' receptor sequences depending on \code{as.list} and \code{group.by}. Each 
#' entry includes CDR3 sequence, V(D)J gene segments, and associated barcodes.
#' 
getIR <- function(input.data,
                  chains,
                  sequence.type = c("aa", "nt"),
                  group.by = NULL,          
                  as.list  = FALSE) {
  
  # Preflight checks-----------------------------------------------------------
  sequence.type <- match.arg(sequence.type)
  col.pos <- if (sequence.type == "aa") "CTaa" else "CTnt"
  
  meta <- if (.is_seurat_or_se_object(input.data)) {
    .grabMeta(input.data)
  } else if (inherits(input.data, "list")) {
    do.call(rbind, input.data)
  } else {
    as.data.frame(input.data)
  }
  rownames(meta) <- meta$barcode %||% rownames(meta)
  
  ok <- c("TRA","TRB","TRG","TRD","Heavy","Light")
  if (anyNA(match(chains, ok)))
    stop("`chains` must be one of: ", paste(ok, collapse = ", "))
  
  if (!is.null(group.by)) {
    if (!group.by %in% names(meta))
      stop("`group.by` column '", group.by, "' not found in metadata.")
    grp_vec <- meta[[group.by]]
  }
  
  # build IR per chain --------------------------------------------------------
  IR <- lapply(chains, function(ch) {
    df <- .get_vdj_matrix(meta[[col.pos]], meta[["CTgene"]], chains)
    df$barcode <- rownames(meta)
    df$chain   <- ch
    if (!is.null(group.by)) df[[group.by]] <- grp_vec
    df
  })
  
  names(IR) <- chains
  
  # assemble return object ----------------------------------------------------
  out <- if (length(chains) == 1L && !as.list) {
    IR[[1L]]
  } else if (!as.list) {
    do.call(rbind, IR)
  } else {             
    IR
  }
  
  if (!is.null(group.by)) {
    # split by grouping variable, preserving list/data.frame structure inside
    split(out, out[[group.by]])
  } else {
    out
  }
}


#' @importFrom SingleCellExperiment colData 
#' @importFrom methods slot
.grabMeta <- function(sc) {
  if (inherits(x=sc, what ="Seurat")) {
    meta <- data.frame(sc[[]], slot(sc, "active.ident"))
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[length(meta)] <- "cluster.active.ident"
    } else {
      colnames(meta)[length(meta)] <- "cluster"
    }
  } else if (inherits(x=sc, what ="SingleCellExperiment")){
    meta <- data.frame(colData(sc))
    rownames(meta) <- sc@colData@rownames
    clu <- which(colnames(meta) == "ident")
    if ("cluster" %in% colnames(meta)) {
      colnames(meta)[clu] <- "cluster.active.idents"
    } else {
      colnames(meta)[clu] <- "cluster"
    }
  }
  return(meta)
}

.split_to_matrix <- function(x, pattern = "_") {
  do.call(rbind, strsplit(x, pattern, fixed = TRUE))
}

#' @keywords internal
.get_vdj_matrix <- function(ct_aa, ct_gene, chain) {
  
  # Pre-process to remove secondary chains (anything after a ';')
  ct_aa_clean <- gsub(";.*", "", ct_aa)
  ct_gene_clean <- gsub(";.*", "", ct_gene)
  
  # Split the strings into alpha/beta chain components.
  aa_list <- strsplit(ct_aa_clean, "_")
  gene_list <- strsplit(ct_gene_clean, "_")
  
  aa_mat <- do.call(rbind, lapply(aa_list, `length<-`, 2))
  gene_mat <- do.call(rbind, lapply(gene_list, `length<-`, 2))
  
  # Determine which column to use based on the requested chain.s.
  is_second_type <- chain %in% c("TRB", "TRD", "Light")
  is_heavy_type <- chain %in% c("TRB", "TRD", "Heavy")
  col_idx <- if (is_second_type) 2 else 1
  
  cdr3_selected <- aa_mat[, col_idx]
  gene_selected <- gene_mat[, col_idx]
  
  # Split the selected gene strings by '.' to get V, D, J, C segments..
  gene_selected_safe <- ifelse(is.na(gene_selected), "", gene_selected)
  gene_parts_list <- strsplit(gene_selected_safe, "\\.")
  
  # Determine the expected number of gene segments and create the matrix.
  max_genes <- if (is_heavy_type) 4 else 3 
  gene_parts_mat <- do.call(rbind, lapply(gene_parts_list, `length<-`, max_genes))
  
  # Assemble the final data frame with the extracted components.
  if (is_heavy_type) {
    df <- data.frame(
      cdr3_aa = cdr3_selected,
      v = gene_parts_mat[, 1],
      d = gene_parts_mat[, 2],
      j = gene_parts_mat[, 3],
      c = gene_parts_mat[, 4],
      stringsAsFactors = FALSE
    )
  } else {
    df <- data.frame(
      cdr3_aa = cdr3_selected,
      v = gene_parts_mat[, 1],
      d = NA, # D gene is not applicable for light/alpha chains
      j = gene_parts_mat[, 2],
      c = gene_parts_mat[, 3],
      stringsAsFactors = FALSE
    )
  }
  
  # Final cleanup
  df[df == "NA" | df == "None" | df == ""] <- NA
  
  return(df)
}
