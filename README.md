# Analysis of *Macrorhabdus ornithogaster* treatment data

This repository contains data and code for the statistical analyses found in the following manuscript:

* Baron, H. R., Leung, C. L., Stevenson, B. C., Gonzalez, M. S., and Phalen, D. N. Evidence of Amphotericin B resistance in *Macrorhabdus ornithogaster* in Australian cage-birds.

All analyses require the free software environment [R](https://www.r-project.org/).

## Chicken treatment trials

R code for the analysis of the chicken treatment trials data can be found in the file `chickens.r`.

### Requirements

This code requires the R package `TMB` and a C++ compiler to be installed.

## Lovebird and budgerigar screening

R code for the analysis of the lovebird and budgerigar screening data can be found in the file `screening.r`.

### Requirements

This code requries the R package `lubridate`.

## Amphotericin B  treatment success 2008–2017

R code for the analysis of the amphotericin B treatment success data can be found in the file `treatment-success.r`.

### Requirements

This code requires the R package `survival`.