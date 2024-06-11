
# xeniumUtils

<!-- badges: start -->
<!-- badges: end -->

This package provides a few utility functions to make handling Xenium In Situ data easier.

The package is currently under development and may be unstable. If there are any problems, bugs, suggestions, or feedback, please feel free to open up an issue.

## Installation

You can install the development version of xeniumUtils from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cindyfang70/xeniumUtils")
```

An example of reading in data:

``` r
library(xeniumUtils)
sfe <- readXenium("Xenium_V1_hPancreas_Cancer_Add_on_FFPE_outs")
```

