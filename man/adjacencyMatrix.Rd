% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjacencyMatrix.R
\name{adjacencyMatrix}
\alias{adjacencyMatrix}
\title{Adjacency matrix from amino acid or nucleotide sequences}
\usage{
adjacencyMatrix(
  input.sequences = NULL,
  normalize = TRUE,
  sequence.dictionary = amino.acids
)
}
\arguments{
\item{input.sequences}{The amino acid or nucleotide sequences to use}

\item{normalize}{Return the values as a function of total number of 
residues (\strong{TRUE}) or frequencies (\strong{FALSE})}

\item{sequence.dictionary}{The letters to use in sequence generation 
(default are all amino acids)}
}
\value{
Adjacency matrix based on input.sequences.
}
\description{
Calculate frequency of adjacency between residues
along a set of biological sequences.
}
\examples{
new.sequences <- generateSequences(prefix.motif = "CAS",
                                   suffix.motif = "YF",
                                   number.of.sequences = 100,
                                   min.length = 8,
                                   max.length = 16)
                          
adj.matrix <- adjacencyMatrix(new.sequences,
                              normalize = TRUE)
                        
}
