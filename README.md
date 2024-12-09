# rsvie

`rsvie` is a simulation package which allows users to evaluate the impact and cost-effectiveness of potential RSV intervention programmes in England and Wales. The user can define intervention program timings, products,  risk of outcomes, costs, and QALYs.

## Installation

First, install [R](https://cran.r-project.org/).

Then, install the development version of 'rsvie' by running the following commands in R.

```
# install.package("devtools") only run if devtools isn't installed.
library(devtools)
github_install("dchodge/rsvie")
library(rsvie)
```

## Quick start and overview of the package

The user must define several components in order to simulate the impact of an RSV intervention programme. These are

1. `add_economics`, the user must define the health outcomes considered, the risk of each health outcome occurring per-infection, and the QALY loss and cost of each outcome.
2. `add_programme`, the user must define the programme characteristics they wish to evaluate, including the timings (start date and end date), the target age group, the product's immunological profile, and the uptake rate.

A detailed vignettes explaining explaining how these work are in the `vignettes` folder.

## Contact

If you have any questions, please email me at david.hodgson@lshtm.ac.uk
