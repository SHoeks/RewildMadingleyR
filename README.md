# MadingleyR Rewilding

MadingleyR version modified to run simulations for the Rewilding Europe project. Two modifications were implemented to the original code in order to improve the realism of the simulations:

-   The competition for resources between herbivore body mass categories was partially lifted by binning herbivores using body mass (into 8 bins) and allowing each bin to feed on independent vegetation stocks.

-   Small-bodied prey (\<150g) were made invisible to predators once the summed cohort density of all cohorts with similar traits dropped below a set threshold (to protect them for going extinct).

The published MadingleyR version can be found [here](https://madingleyr.github.io/MadingleyR/). Modifications were made to MadingleyR version v1.0.0 and source code version v2.00 published on the 11th of February 2022.

### Installation

Installing the package into your R library can be done using the following code. Please note, that the source code is only compiled for Linux. Instructions for how to compile the code for Windows and Mac can be found [here](https://madingleyr.github.io/MadingleyR/SourceCode/CompileMac.html) and [here](https://madingleyr.github.io/MadingleyR/SourceCode/CompileWindows.html).

``` r
# Download the zip from this repository
# link: 

# Install package from the local zip using devtools
devtools::install_local("MadingleyRewilding-main.zip")
```

### Usage

Usage of this version of the package is identical to the published MadingleyR package with the exception of three input arguments:

-   XXX1: this controls the body mass below which the hetetrophic cohorts can become invisible once their density drops below the threshold set by XXX2.

-   XXX2: this defines the minimum density threshold (in individuals/km^2^) below which the small-bodied heterotrophs (with a body mass \< XXX1) will become invisible.

-   XXX3: a vector defining the borders of the body masses bins for herbivores (in grams). Herbivores between these borders will be binned and thus competing for the same resources. The number of bin borders defined here will determine the number of herbivore groups (*n_bin_borders* - 1 = *n_herbivore_groups*), the number of herbivore groups will be used to divide the vegetation over equal pools.

Below a snippet of code to illustrate the use of the three input arguments described above.

``` r
# Source the package
library(MadingleyRewilding)
```
