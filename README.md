<!-- badges: start -->
[![Build Status](https://travis-ci.org/mrc-ide/vivax.svg?branch=master)](https://travis-ci.org/mrc-ide/vivax)[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/mrc-ide/vivax?branch=master&svg=true)](https://ci.appveyor.com/project/mrc-ide/vivax)[![codecov](https://codecov.io/github/mrc-ide/vivax/branch/master/graphs/badge.svg)](https://codecov.io/github/mrc-ide/vivax)
<!-- badges: end -->

# P Vivax individual based model

This is an R wrapper for Dr. Michael White's model for P. Vivax malaria
transmission.

This model was first described in the article "Impact of expanding levels of
malaria control interventions on Plasmodium vivax: a mathematical modelling study".

## Installation

The package can be installed from github using the "remotes" library

```R
library('remotes')
install_github('mrc-ide/vivax')
```

For development it is most convenient to run the code from source. You can
install the dependencies in RStudio by opening the project and selecting "Build" > "Install and Restart"

Command line users can execute:

```R
library('remotes')
install_deps('.', dependencies = TRUE)
```

## Usage

To run the simulation with the default parameters you
can execute the following code:

```R
library('vivax')
output <- run_simulation()
```
## Parameterisation

You can find the default paramters
[here](https://github.com/mrc-ide/vivax/tree/master/inst/defaults).

To override either the model, farauti, koliensis, punctulatis parameters by
passing named lists to the `run_simulation` function as below:

```R
output <- run_simulation(model=list(EIR_equil=.1), koliensis=list(mu_P=.5, dry_seas=.2))
```

### Interventions

The available interventions are documented
[here](https://github.com/mrc-ide/vivax/tree/master/R/interventions.R).

To override interventions pass your named list of parameters to `run_simulation` function as below:

```R
output <- run_simulation(interventions = list(LLIN_years=c(1995, 1996), LLIN_cover=c(.6, .8)))
```

## Code organisation

*src/model.cpp* - this is the original model code. `run_simulation_from_path` is the
entry point from R.

*R/run.R* - contains code for parameterisation and running of the model.

*tests* - are divided into unit and integration tests.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to
discuss what you would like to change.

Please make sure to update tests and documentation as appropriate.

Code reviews will be carried out in-line with RESIDE-IC's [review
process](https://reside-ic.github.io/articles/pull-requests/)

## License
[MIT](https://choosealicense.com/licenses/mit/)
