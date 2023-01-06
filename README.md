# Analysis of *Macrorhabdus ornithogaster* treatment data

This repository contains data and code for the statistical analyses found in the following manuscripts:

* Baron, H. R., Stevenson, B. C., and Phalen, D. N. (2021) Comparison of in-clinic diagnostic testing methods for *Macrorhabdus ornithogaster*. *Journal of Avian Medicine and Surgery*, *35*(1), 37&ndash;44.

* Baron, H. R., Stevenson, B. C., and Phalen, D. N. (2020) Inconsistent efficacy of water soluble Amphotericin B for the treatment of *Macrorhabdus ornithogaster* in a budgerigar (*Melopsittacus undulatus*) aviary. *Australian Veterinary Journal*, *98*(7), 333&ndash;337.

* Baron, H. R., Leung, C. L., Stevenson, B. C., Gonzalez, M. S., and Phalen, D. N. (2019) Evidence of Amphotericin B resistance in *Macrorhabdus ornithogaster* in Australian cage-birds. *Medical Mycology*, *57*(4), 421&ndash;428.

All analyses require the free software environment [R](https://www.r-project.org/).

## Comparison of in-clinic diagnostic testing methods for *Macrorhabdus ornithogaster*

R code for the analysis of the point-of-care diagnostic modality comparison can be found in the file `clinic-modality-comparison/comparison.r`.

Requirements:
* The R package `TMB`
* A C++ compiler.

## Inconsistent efficacy of water soluble Amphotericin B for the treatment of *Macrorhabdus ornithogaster* in a budgerigar (*Melopsittacus undulatus*) aviary

R code for the analysis of the budgerigar faecal shedding data can be found in the file `treatment-effects/treatment.r`.

Requirements:
* The R package `TMB`
* A C++ compiler

## Evidence of Amphotericin B resistance in *Macrorhabdus ornithogaster* in Australian cage-birds

### Chicken treatment trials

R code for the analysis of the chicken treatment trials data can be found in the file `ab-resistance/chickens.r`.

Requirements:
* The R package `TMB`
* A C++ compiler

### Lovebird and budgerigar screening

R code for the analysis of the lovebird and budgerigar screening data can be found in the file `ab-resistance/screening.r`.

Requirements:
* The R package `lubridate`

### Amphotericin B treatment success 2008–2017

R code for the analysis of the Amphotericin B treatment success data can be found in the file `ab-resistance/treatment-success.r`.

Requirements:
* The R package `survival`
