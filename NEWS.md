# immApex VERSION 1.7.2

## BUG FIXES
* Fixed `sequenceEncoder()` and `sequenceDecoder()` property mode silently mismatching amino-acid identity to property values for `crucianiProperties`, `MSWHIM`, and `ProtFP`. These property sets' column order in `Peptides::AAdata` is not already canonical (MSWHIM/ProtFP are stored alphabetically; crucianiProperties has an isolated E/Q swap), and `.aa.property.matrix()` was returning that raw storage order instead of reordering to match the canonical amino-acid order its callers assume positionally. For example, encoding "R" under MSWHIM silently returned Cysteine's values. `.aa.property.matrix()` now always returns canonical-AA-ordered columns, matching `calculateProperty()`'s existing contract and the documented `sequenceEncoder()`/`sequenceDecoder()` interface. Verified against `Peptides`' own independent public scoring functions (`mswhimScores()`, `crucianiProperties()`, `protFP()`).

# immApex VERSION 1.7.1

## PERFORMANCE
`buildNetwork()` now scales to large, clonally expanded repertoires. The C++ engine was rewritten around candidate generation instead of all-pairs comparison, with no new dependencies and identical edge output for the default settings:

* Identical sequences collapse to one representative before any distance is computed, so clonal expansion no longer drives quadratic work.
* Length and V/J blocking act as an index. Only length-compatible representatives in the same block are compared.
* Levenshtein and Damerau use deletion-neighborhood (SymSpell) hashing for small edit distances, with Damerau indexed at radius `2k`. Hamming uses pigeonhole segment indexing. Both fall back to length-blocked comparison outside their range.
* Edges accumulate as integer node indices and gain their labels in R, which avoids building millions of barcode strings inside the parallel core.

Speedups range from roughly 2x on mixed repertoires with many V genes to about 70x when one block holds many unique sequences at a tight threshold.

## NEW FEATURES
* Added the `expand` argument to `buildNetwork()`. `"clique"` (default) materializes every pairwise edge and reproduces the exact edge multiplicity that community-detection clustering expects. `"star"` links identical sequences through a single hub and connects related groups hub to hub. `"star"` produces far fewer edges and preserves connected components exactly, so it is a large memory win when the downstream step depends on connectivity.

## BUG FIXES
* Fixed a banded Levenshtein error that dropped pairs whose only optimal alignment rides the band edge, including identical sequences at a tight normalized threshold. The boundary cell read the wrong diagonal value. The fix only recovers true edges and never removes valid ones.
* Made the normalized-threshold cutoff robust to floating-point rounding so a normalized distance exactly equal to the threshold is kept, matching the documented inclusive behavior.

 # immApex VERSION 1.5.4

## BUG FIXES
* Fixed vignette build failure when IMGT network data is unavailable (e.g., on Bioconductor build servers)
* Fixed operator precedence bug in `inferCDR()` reference validation check

# immApex VERSION 1.5.3

## UNDERLYING CHANGES
* `getIMGT()` now uses [immReferent](https://github.com/BorchLab/immReferent) as a backend for downloading and caching IMGT reference sequences
* Removed `frame`, `max.retries`, and `verbose` parameters from `getIMGT()`
* Added `refresh` parameter to `getIMGT()` for controlling cache behavior
* Removed `hash`, `httr`, and `rvest` package dependencies (now handled by immReferent)
* Replaced `.parseSpecies()` with `.mapSpecies()` for immReferent species name compatibility

# immApex VERSION 1.5.2

- improved error message for the `propertyEncoder()` function when the suggested `Peptides` package is not installed.

# immApex VERSION 1.5.0

- Bioconductor devel branch update to match versions

# immApex VERSION 1.4.3

## NEW FEATURES
The `buildNetwork()` function has been significantly enhanced to include:

* **Multiple distance metrics**: Added support for five distance metrics via the `dist_type` parameter:
  - `"levenshtein"` (default): Standard edit distance
  - `"hamming"`: Substitutions only (requires equal-length sequences)
  - `"damerau"`: Levenshtein with transpositions
  - `"nw"`: Needleman-Wunsch global alignment
  - `"sw"`: Smith-Waterman local alignment

* **Flexible normalization**: New `normalize` parameter supports three modes:
  - `"none"` (default): Raw distance values (backward compatible)
  - `"maxlen"`: Normalize by max(length(seq1), length(seq2))
  - `"length"`: Normalize by mean sequence length

# immApex VERSION 1.3.7

## UNDERLYING CHANGES
* Ensure `buildNetwork()` returns a symmetric distance matrix and drops nonstructural zeros

# immApex VERSION 1.3.6

## UNDERLYING CHANGES
* Rebuilt ```propertyEncoder()``` and ```onehotEncoder()``` to use C++ backend with encodeSequences.cpp
* Rebuilt ```buildNetwork()``` via C++ integration (fastEditEdges.cpp)
* Motif quantification in ```calculateMotif()``` use C++ (calculateMotif.cpp)
* Reduce overall dependencies on external packages
* Improved speed of ```generateSequences()```, ```mutateSequences()```, ```adjacencyMatrix()```, ```tokenizeSequences()```, ```geometricEncoder```
* Added group.by argument to ```getIR()```
* Defunct ```variationalSequences()``` 
* Removed keras3/tensorflow suggests/dependencies
* Removed baslisk environment, python no longer required
* Updated C++ backend to C++17
* Robust checks before running `getIMGT()` for website functioning
* Removing special unicode symbols

## NEW FEATURES
* ```calculateEntropy()``` added to calculate the positional entropy along a biological sequence
* Added basic diversity metrics (exported) to support ```calculateEntropy()```
* ```calculateFrequency()``` added to calculate the positional frequency along a biological sequence
* ```calculateMotif()``` added to get motif quantification of sequences
* ```calculateGeneUsage()``` added for single/paired gene enumeration
* Added ```scaleMatrix()``` for comprehensive scale/transformation functions
* Added ```summaryMatrix()``` for fast summarization of matrix values

## BUGS
* Fixed NT position issue with ```inferCDR()```
* Fixed ```generateSequence()``` range issue for single-length sequences
* Fixed binary frequency issue in ```positionalEncoder()```
* Fixed ```buildNetwork()``` for relative thresholding returns relative value

# immApex VERSION 1.2.2
* Updated unit tests and vignette check
* Converted package to basilisk from reticulate

# immApex VERSION 1.2.0
* Update version for BioConductor release

# immApex VERSION 1.0.5

## UNDERLYING CHANGES
* ```getIMGT()``` checks for availability of IMGT website
* Expanded Unit Tests
* Unit Tests and Vignette now evaluate for proper python installation overall

# immApex VERSION 1.0.4

## UNDERLYING CHANGES
* Optional testthat ```variationalSequences()``` evaluate presence of Keras

# immApex VERSION 1.0.3

## UNDERLYING CHANGES
* Drop evaluation of ```variationalSequences()``` example

# immApex VERSION 1.0.2

## UNDERLYING CHANGES
* Vignette includes eval of keras installation for certain chunks

# immApex VERSION 1.0.1

## UNDERLYING CHANGES
* Fix issue with optimizer call in ```variationalSequences()```
* Move package to keras3
* Version Bump to be consistent with Bioconductor release

# immApex VERSION 0.99.5

## UNDERLYING CHANGES
* Remove keras installation in vignette
* Removed lazydata loading
* Added @return to the data manual entries
* Update stop message for ```propertyEncoder()```
* Removed set.seed from ```variationalSequences()```
* Added runnable example for ```variationalSequences()```
* Added better handling of no connection in ```getIMGT()```
* Added **max.retries** to ```getIMGT()```
* Updated testthat for new IMGT version

# immApex VERSION 0.99.3

## UNDERLYING CHANGES
* Defining python version 3.10 for keras installation

# immApex VERSION 0.99.2

## UNDERLYING CHANGES
* Added discrete check and install_keras() to vignette


# immApex VERSION 0.99.1

## UNDERLYING CHANGES
* Added warning to ```getIMGT()``` to mention IMGT license on first use
* import magrittr instead of dplyr
* Fixed issue with license error
* Added example to ```positionalEncoder()```
* Replaced discrete calling of "es =" in vignette
* Update readme to include keras based installation
* added set.seed() to vignette
* added introduction section to vignette with mention of other packages
* removed installation section from vignette (can't add BioC instructions until accepted)
* Fixed rendering issue with the vignette
* Subsampled example data and removed Omniscope example
* Updated testthat for ```formatGenes()``` and ```inferCDR()``` for new example data
* Reduced size of test data
* Minimized usage of sapply() and X:Y per BiocCheck suggestions
* Modified print() messages to messages()
* Add verbose parameter to turn on/off processing steps
* Added links for example data formats
