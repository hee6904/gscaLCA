## Test environments
* local OS X install, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* We received feedback on our R package version 0.0.1. After revise, we resubmit it as version 0.0.2

* Win-builder generated 1 NOTE: * possibly misspelled words in DESCRIPTION, but the spelling is correct.

* The 'doSNOW' package cannot be replaced with other packages in order to show the progress bar with multiple core running. 

* The examples in this package use a single core. According to preference by users, multiple cores can be used. 

* This packcage provides both an object of information and printed output in console.
We think that printing out the output in the console is crucial, and it will be convenient for users.
