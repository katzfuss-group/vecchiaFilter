# VEnKF

This package is based on my research project with Dr. Matthias Katzfuss using Vecchia approximation to regularize covariance estimation within an Ensemble Kalman Filter update step. The package will contain hidden functions that perform background Vecchia calculations based upon discussions with Brian Kidd. The main purpose of this package will instead be focused upon the update function based on implementing all of these Vecchia functions to calculate the EnKF update step. The package will also contain several other update functions to compare these update steps with example graphical representations of possible simulations. Finally, there will be an implementation of a Lorenz Chaos Model to run numerical simulations of the update functions. The package also contains several small example datasets. 

## Installation with vignettes

```{r eval = F}
devtools::install_github("wboyles/VEnKF", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
library(VEnKF)
```


