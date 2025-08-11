# immApex VERSION 1.3.7

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
