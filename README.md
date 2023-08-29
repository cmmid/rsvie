# rsvie

`rsvie` is a simulation package which allows users to evaluate the impact and cost-effectiveness of potential RSV intervnetion programmes in England and Wales. The user can define interventin nprogramme timings, products, change the risk of outcomes, costs, and QALYs.

## Installation

First install [R](https://cran.r-project.org/).

Then install the development version of 'rsvie' by running the following commands in R.

```
# install.package("devtools") only run if devtools isn't installed.
library(devtools)
github_install("dchodge/rsvie")
library(rsvie)
```

## Quick start and overview of package

The user must define several components in order to simulate the impact of an RSV intervention programme. These are
1. `add_economics`, the user must define the health outcomes considered, the risk of each healthcome occuring per-ifection, and the QALY loss and cost of each outcome.
2. `add_programme`, the user must define the characteristics of the programme they wish the evaluate, including the timings (strat date and end date), the target age-group, the immunological profile of the product, and the rate of uptake. 

A detailed vignette explaning how to define these characteristics can be found in `vignettes/example.Rmd`. 

## Example vignettes

We also include several vignettes showing how to evaluate various programmes. These are found in `vignettes/`.

## Contact

Any questions, please email me on david.hodgson@lshtm.ac.uk