# hmm-missing-data-paper

Code and RMarkdown for "Ignorable and non-ignorable missing data in hidden Markov models"

## Packages and versions

For R version and package versions used to compile the paper, please see the `sessionInfo.txt` file.

## Paper

Text and code for the paper is contained in the `hmm-missing-data-paper.Rmd` R Markdown file.
This file contains all the code for the analyses, apart from the simulations. 
Code for the simulations is in separate files (see below).

Note that the hidden Markov analyses in chunks `fit-hmms` and 
especially `constrained-state-dependent` require substantial
computation time.

The code in the `load-data` chunk retrieves the data from the Schizophrenia
study from https://hedeker.people.uic.edu/SCHIZREP.DAT.txt. If that dataset 
is no longer retrievable from this address, a copy should be available from
https://web.archive.org/web/20100731195849/http://tigger.uic.edu/~hedeker/SCHIZREP.DAT.txt

## Simulations

Code for running the simulation studies is provided in the files `simulation1.R` 
to `simulation5.R`. Also provided are the results, in the files `simulation1.RData`
to `simulation5.RData`, as well as console output whilst running the simulations, in 
files `simulation.Rout` to `simulation5.Rout`. These latter files also contain the
output of `proc.time()` which shows the duration.

Simulations to compute the maximum classification accuracy for the HMMs in provided in the `classification_accuracy.R` file. The results are provided in the `classification_accuracy.RData` file.

## Reproducing the results

* Source the files `simulation1.R`, `simulation2.R`, `simulation4.R`, `simulation5.R`
* Source the file `classification_accuracy.R`
* Run e.g. `knitr::knit("hmm-missing-data.paper.Rmd`)

