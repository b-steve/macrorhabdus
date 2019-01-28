# Analysis of *Macrorhabdus ornithogaster* treatment data

This repository contains data and code for the statistical analyses found in the following manuscripts:

* Baron, H. R., Leung, C. L., Stevenson, B. C., Gonzalez, M. S., and Phalen, D. N. (in press) Evidence of Amphotericin B resistance in *Macrorhabdus ornithogaster* in Australian cage-birds. *Medical Mycology*.

* Baron, H. R., Stevenson, B. C., and Phalen, D. N. (in preparation) Comparison of point-of-care diagnostic modalities for *Macrorhabdus ornithogaster*.

All analyses require the free software environment [R](https://www.r-project.org/).

## Chicken treatment trials

R code for the analysis of the chicken treatment trials data in Baron et al. (in press) can be found in the file `ab-resistance/chickens.r`.

### Requirements

This code requires the R package `TMB` and a C++ compiler.

## Lovebird and budgerigar screening

R code for the analysis of the lovebird and budgerigar screening data in Baron et al. (in press) can be found in the file `ab-resistance/screening.r`.

### Requirements

This code requries the R package `lubridate`.

## Amphotericin B  treatment success 2008–2017

R code for the analysis of the amphotericin B treatment success data in Baron et al. (in press) can be found in the file `ab-resistance/treatment-success.r`.

### Requirements

This code requires the R package `survival`.

## Comparison of point-of-care diagnostic modalities

R code for the analysis of the point-of-care diagnostic modality comprison in Baron, Stevenson, and Phalen (in preparation) can be found in the file `modality-comparison/comparison.r`.

### Requirements

This code requires the R packates `TMB` and a C++ compiler.