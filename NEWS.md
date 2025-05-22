# immApex VERSION 1.2.1
* Update unit testing 
* Removed src .o and .so files

# immApex VERSION 1.0.6
* Update unit testing and keras vignette check

# immApex VERSION 1.0.5

## UNDERLYING CHANGES
* ```getIMGT()``` checks for availability of IMGT website
* Expanded Unit Tests
* Unit Tests and Vignette now evaluate for proper python installation overall

# immApex VERSION 1.0.4

## UNDERLYING CHANGES
* Optional testthat variationalSequences() evaulate presence of Keras

# immApex VERSION 1.0.3

## UNDERLYING CHANGES
* Drop evaluation of variationalSequences() example

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
