% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probabilityMatrix.R
\name{probabilityMatrix}
\alias{probabilityMatrix}
\title{Position Probability Matrix for Amino Acid or Nucleotide Sequences}
\usage{
probabilityMatrix(
  input.sequences,
  max.length = NULL,
  convert.PWM = FALSE,
  background.frequencies = NULL,
  sequence.dictionary = amino.acids,
  padding.symbol = ".",
  verbose = TRUE
)
}
\arguments{
\item{input.sequences}{The amino acid or nucleotide sequences to use}

\item{max.length}{Additional length to pad, NULL will pad sequences 
to the max length of input.sequences}

\item{convert.PWM}{Convert the matrix into a positional weight matrix 
using log likelihood}

\item{background.frequencies}{Provide amino acid or nucleotide frequencies
for the positional weight matrix. If NULL, assumes uniform likelihood.}

\item{sequence.dictionary}{The letters to use in sequence generation 
(default are all amino acids)}

\item{padding.symbol}{Symbol to use for padding at the end of sequences}

\item{verbose}{Print messages corresponding to the processing step}
}
\value{
A matrix with position specific probabilities or weights
}
\description{
Use this to generate a position-probability or weight matrix 
for a set of given sequences.
}
\examples{
new.sequences <- generateSequences(prefix.motif = "CAS",
                                   suffix.motif = "YF",
                                   number.of.sequences = 100,
                                   min.length = 8,
                                   max.length = 16)
                          
PPM.matrix <- probabilityMatrix(new.sequences)
                        
}
