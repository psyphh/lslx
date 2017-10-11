# lslx
lslx is a package for fitting semi-confirmatory structural equation modeling (SEM) via penalized likelihood (PL) developed by Huang, Chen, and Weng (2017). In this semi-confirmatory method, an SEM model is distinguished into two parts: a comfirmatory part and an exploratory part. The confirmatory part includes all of the freely estimated parameters and fixed parameters that are allowed for theory testing. The exploratory part is composed by a set of penalized parameters describing relationships that cannot be clearly determined by available substantive theory. By implementing a sparsity-inducing penalty and choosing an optimal penalty level, the relaitonships in the exploratory part can be efficiently determined by the sparsity pattern of these penalized parameters. 

# Installation
The following code will install `lslx` and all the dependencies. 
``` r
devtools::install_github("psyphh/lslx", 
                          force = TRUE, 
                          dependencies = TRUE,
                          build_vignettes = TRUE)
```
If you have not install package `devtools` before, please install `devtools` first.
For windows users, `Rtools` (https://cran.r-project.org/bin/windows/Rtools/) should be installed first.

# Usage
`lslx` is an R6ClassGenerator for constructing an `lslx` object that has methods for fitting semi-confirmatory SEM. In a simpliest case, the use of `lslx` involves three major steps
1. Initialize a new `lslx` object by specifying a model and importing a data set.
``` r
r6_lslx <- lslx$new(model, data)
```
2. Fit the specified model to the imported data with specified fitting control.
``` r
r6_lslx$fit(penalty_method, lambda_grid, gamma_grid)
```
3. Summarize the fitting results with specified selector.
``` r
r6_lslx$summarize(selector)
```

# Tutorial
You can learn how to use `lslx` by examples with the vignettes.
``` r
vignette("lslx-example-1")
vignette("lslx-example-2")
vignette("lslx-example-3")
vignette("lslx-example-4")
```
Of course, help file is also available.
``` r
?lslx
```

# On Going Jobs
The current beta version only provides limited functionality. 
The road map for major updating is
1. Implementing valid post model selection inference methods.
2. Implementing pairwise likelihood approach for ordinal data.

# Bug Reporting
If you find any bug, please mail psyphh@gmail.com.
