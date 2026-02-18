
# Install dependencies

In order to use `espatsmo` to its fullest extent, users should install INLA. This will bring the possibility of fitting a log-Gaussian Cox point process model as implemented in the ppmLGCP function. However, INLA has to be installed manually using the following commands:

`install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)`

Then, linux users should load INLA, and then:

`inla.binary.install()`

# Install `espatsmo` 

To install `espatsmo`, use:

`devtools::install_github("gerardommc/espatsmo")`