# SAZED

## Overview
SAZED (Spectral and Average Autocorrelation Zero Distance Density) 
is a package and ensemble method for estimating the season length of 
a seasonal time series. SAZED is aimed at practitioners, as SAZED 
employs only domain-agnostic pre-processing and does not depend on 
parameter tuning or empirical constants.

## Installation
This package requires the debian packages libxml2-dev, libfftw3-dev 
and zlib1g-dev. Once you have those, then simply execute:

``` r
install.packages("sazed")
```


## Usage
Estimate the season length of a seasonal time series with the 
ensemble method SAZED:

``` r
library(sazed)

season_length <- 26
y <- sin(1:400*2*pi/season_length)
sazed(y)
```

All components of the SAZED ensemble are also available separately. 
For more information on them, see the package manual and examples 
therein.
