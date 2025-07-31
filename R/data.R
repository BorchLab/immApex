#' A list of IMGT gene names by genes, loci, and species
#' 
#' @description A list of regularized gene nomenclature to use for 
#' converting for data for uniformity. Data is organize by gene region, 
#' loci and species. Not all species are represented in the data and
#' pseudogenes have not been removed.
#' 
#'  
#' @docType data
#' @usage data("immapex_gene.list")
#' @name immapex_gene.list
#' @return List of gene nomenclature by region, loci, and species. 
#' 
NULL

#' Example contig data for Apex
#' 
#' @description Contains a collection of bulk or paired TCR sequences in the respective formats in the form of a 
#' list from the following sources:
#' 
#' \itemize{
#'   \item{TenX: 10k_Human_DTC_Melanoma_5p_nextgem_Multiplex from \href{https://www.10xgenomics.com/datasets/10k-human-dtc-melanoma-NextGEM-5p}{10x Website}.}
#'   \item{AIRR: Human_colon_16S8157851 from \href{https://pubmed.ncbi.nlm.nih.gov/37055623/}{PMID: 37055623}.}
#'   \item{Adaptive: Adaptive_2283_D0 from \href{https://pubmed.ncbi.nlm.nih.gov/36220826/}{PMID: 36220826}.}
#' }
#'
#' More information on the data formats are available: \href{https://docs.airr-community.org/en/stable/}{AIRR}, 
#' \href{https://clients.adaptivebiotech.com/assets/downloads/immunoSEQ_AnalyzerManual.pdf}{Adaptive}, and 
#' \href{https://www.10xgenomics.com/support/single-cell-immune-profiling}{TenX}.
#' 
#' @docType data
#' @usage data("immapex_example.data")
#' @name immapex_example.data
#' @return List of 3 example data sets for 10x, AIRR and Adaptive contigs.
#' 
NULL

#' List of amino acid substitution matrices
#' 
#' @description A list of amino acid substitution matrices, using the Point 
#' Accepted Matrix (PAM) and BLOck SUbstitution Matrix (BLOSUM) approaches. 
#' A discussion and comparison of these matrices are available at
#' \href{https://pubmed.ncbi.nlm.nih.gov/21356840/}{PMID: 21356840}.
#' 
#' \itemize{
#'  \item{BLOSUM45}
#'  \item{BLOSUM50}  
#'  \item{BLOSUM62}
#'  \item{BLOSUM80}
#'  \item{BLOSUM100}
#'  \item{PAM30}
#'  \item{PAM40}
#'  \item{PAM70}    
#'  \item{PAM120}
#'  \item{PAM250}
#'  }   
#' 
#' @docType data
#' @name immapex_blosum.pam.matrices
#' @usage data("immapex_blosum.pam.matrices")
#' @return List of 10 substitution matrices
#' 
NULL