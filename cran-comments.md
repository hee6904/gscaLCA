## Test environments
* local OS X install, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is the second submission after revise.
* Win-builder generated 1 NOTE: * possibly misspelled words in DESCRIPTION, but the spelling is correct

* The 'doSNOW' package cannot substituted with other packages in order to show progress bar with multiple core running. 

* The examples in this package uses a single core. According to preference by users, multiple core can be used. 

* Our program provides information as both forms of an object and printed output in console. 
We think that printing out the output in the console is really important. It will be convenient for users.
