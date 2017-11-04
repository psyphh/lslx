## What is lslx?
**lslx** is a package for fitting semi-confirmatory structural equation modeling (SEM) via penalized likelihood (PL) with lasso or minimax concave penalty (MCP) developed by Huang, Chen, and Weng (2017) <<[doi:10.1007/s11336-017-9566-9](doi:10.1007/s11336-017-9566-9)>>. In this semi-confirmatory method, an SEM model is distinguished into two parts: a confirmatory part and an exploratory part. The confirmatory part includes all of the freely estimated parameters and fixed parameters that are allowed for theory testing. The exploratory part is composed by a set of penalized parameters describing relationships that cannot be clearly determined by available substantive theory. By implementing a sparsity-inducing penalty and choosing an optimal penalty level, the relationships in the exploratory part can be efficiently determined by the sparsity pattern of these penalized parameters. 

**lslx** can be also seen as a package for conducting usual SEM with several state-of-art inference methods, including sandwich standard error formula, mean-adjusted likelihood ratio test, and two step method with auxiliary variables for missing data. **lslx** also supports multi-group analysis for evaluating group heterogeneity. For now, the major limitations of **lslx** are that (1) it cannot impose linear or non-linear constraints for coefficients; (2) it cannot make valid inference under clustered or dependent data; and (3) it cannot handle the presence of ordinal responses.

## Usage
`lslx` is an R6ClassGenerator for constructing an `lslx` object that has methods for fitting semi-confirmatory SEM. In a simplest case, the use of `lslx` involves three major steps

1. Initialize a new `lslx` object by specifying a model and importing a data set.
``` r
r6_lslx <- lslx$new(model, data)
```
2. Fit the specified model to the imported data with specified fitting control.
``` r
r6_lslx$fit(penalty_method, lambda_grid, delta_grid)
```
3. Summarize the fitting results with specified selector.
``` r
r6_lslx$summarize(selector)
```

## Installation
To install the [CRAN](https://CRAN.R-project.org/package=lslx) version, just simply type
``` r
install.packages("lslx", dependencies = TRUE)
```

## Bug Reporting
If you find any bug, please post it on [issue](https://github.com/psyphh/lslx/issues).