# Analysis of *Macrorhabdus ornithogaster* treatment data

This repository contains data and code for the statistical analyses found in the following two manuscripts:

* Baron, H. R., Leung, C. L., Stevenson, B. C., Gonzalez, M. S., and Phalen, D. N. Treatment trials with amphotericin B for the treatment of *Macrorhabdus ornithogaster in experimentally infected white leghorn chickens (*Gallus gallus domesticus*).

* Baron, H. R., Leung, C. L., Stevenson, B. C., Gonzalez, M. S., and Phalen, D. N. Long-term monitoring of faecal shedding of *Macrorhabdus ornithogaster* following treatment with amphotericin B and analysis of treatment success with amphotericin B against *Macrorhabdus ornithogaster* 2008–2017.

All analyses require the free software environment [R](https://www.r-project.org/).

## Chicken treatment trials

R code for the analysis of the chicken treatment trials data (described in the first manuscript) can be found in the file `chickens.r`.

### Requirements

This code requires the R package `TMB` and a C++ compiler to be installed.

## Lovebird and budgerigar screening

R code for the analysis of the lovebird and budgerigar screening data (described in the second manuscript) can be found in the file `screening.r`.

## Amphotericin B  treatment success 2008–2017

R code for the analysis of the amphotericin B treatment success data (described in the second manuscript) can be found in the file `treatment-success.r`.