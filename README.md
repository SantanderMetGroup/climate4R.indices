## What is `climate4R.indices`?

**climate4R** is a set of R packages for transparent climate data access, post processing (including bias correction and downscaling) and visualization. For more information and references, visit the [climate4R page](http://www.meteo.unican.es/climate4r).

`climate4R.indices` is the package to compute several indices within the climate4R framework, therefore is seamlessly integrated with the **climate4R** data structures, and provides support for parallel computing.


****

### Installation

The recommended procedure for installing the package is using the devtools package. Note that this package depends on [`transformeR`](https://github.com/SantanderMetGroup/transformeR), another package from the **climate4R** bundle. Thus:

```R
devtools::install_github(c("SantanderMetGroup/transformeR", "SantanderMetGroup/climate4R.indices"))
```

A list of all available impact indices and the atomic functions calculating them is printed on screen with:

```R
library(climate4R.indices)
indexShow()
?indexGrid   # see the examples 
```

A set of circulation indices is also available:
```R
library(climate4R.indices)
circIndexShow()
?circIndexGrid   # see the examples 
```

Reference and further information: 

**[General description of the climate4R framework]** Iturbide et al. (2019) The R-based climate4R open framework for reproducible climate data access and post-processing. **Environmental Modelling and Software**, 111, 42-54. https://doi.org/10.1016/j.envsoft.2018.09.009
Check out the companion notebooks for the two examples [GitHub](https://github.com/SantanderMetGroup/notebooks).

