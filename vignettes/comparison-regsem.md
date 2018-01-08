Comparison with regsem Package
================

[**regsem**](https://CRAN.R-project.org/package=regsem) is a package for fitting regularized structural equation modeling. It can also implements lasso penalty to obtain sparse estimate.

In this document, we try to compare the fitting results made by **lslx** and **regsem**. The example in **regsem** is adopted, which fits the data of Holzinger and Swineford (1939) by a one-factor model. Our experiment shows that

-   **lslx** and **regsem** yield slightly different results;
-   **lslx** is about 20 times faster than **regsem**.

Since the results of **lslx** and **regsem** are inconsistent, [**lsl**](https://CRAN.R-project.org/package=lsl) is also used to fit the data. We found that **lslx** and **lsl** yield numerically the same results. By the fact that the objective function values made by **lslx** and **lsl** are smaller then that made by **regsem**, we think that the solution made by current version of **regsem** (0.9.2) can be improved.

Note that our comparison is made based on the following versions of packages. The comparison result can be different if other versions are used.

``` r
packageVersion("lslx"); packageVersion("regsem"); packageVersion("lsl")
```

    [1] '0.6.0.9002'

    [1] '1.0.6'

    [1] '0.5.6'

### Syntax for Analysis

We first compare the code syntax of **lslx** and **regsem** by using the example in **regsem**. The code of **lslx** is

``` r
library(lslx)
model_lslx <-
"f1 :=> fix(1) * x1 + pen() * x2 + pen() * x3 + x4 + x5 + x6 + pen() * x7 + pen() * x8 + pen() * x9"

r6_lslx <- lslx$new(model = model_lslx,
                    data = lavaan::HolzingerSwineford1939,
                    verbose = FALSE)

r6_lslx$fit_lasso(lambda_grid = 0.2, verbose = FALSE)
```

In this code, a one-factor model is specified with five loadings to be penalized. **lslx** fits the model to data with lasso penalty with `lambda = 0.2` - only single penalty level is considered for simplicity.

Now we present the code of **regsem** as follows

``` r
library(regsem)
```

    Warning: package 'regsem' was built under R version 3.4.3

    Warning: package 'Rcpp' was built under R version 3.4.3

``` r
model_lavaan <- 
"f =~ 1 * x1 + l1 * x2 + l2 * x3 + l3 * x4 + l4 * x5 + l5 * x6 + l6 * x7 + l7 * x8 + l8 * x9"

fit_lavaan <- lavaan::cfa(model = model_lavaan, 
                          data = lavaan::HolzingerSwineford1939, 
                          meanstructure = TRUE)

fit_regsem <- regsem(fit_lavaan, lambda = 0.2, type = "lasso",
                     optMethod = "rsolnp",
                     pars_pen = c("l1", "l2", "l6", "l7", "l8"),
                     solver = TRUE, quasi = TRUE,
                     solver.maxit = 100, line.search = TRUE)
```

The `regsem()` function uses `lavaan` object as input for further regularized fitting. Hence, model specification in **regsem** relies on **lavaan** syntax, that is attractive for **lavaan** users. `lambda = 0.1` in **regsem** is theoretically equivalent to `lambda = 0.2` in **lslx** (see the discussion in section of Fitting Results). `type = "lasso"` specifies the penalty method and `pars_pen = c("l1", "l2", "l6", "l7", "l8")` specifies which coefficients should be penalized. Other arguments are used for coordinate descent.

### Fitting Results

The objective function values made by **lslx** and **regsem** are

``` r
# objective function value in lslx (on lslx scale)
r6_lslx$extract_numerical_condition()["objective_value"]
```

    objective_value 
           1.246158 

``` r
# objective function value in regsem (on regsem scale)
fit_regsem$optim_fit
```

    [1] 0.6247178

The function values made by **lslx** and **regsem** seem to be on different scale. We already know the scale of objective function in **lslx**. To understood the interpretation of the objective function value in **regsem**, we first reconduct analysis with `lambda = 0` and see the objective function value and loss function value

``` r
fit_regsem_0 <- regsem(fit_lavaan, lambda = 0, type = "lasso",
                       pars_pen = c("l1", "l2", "l6", "l7", "l8"))
# objective function value in regsem (on regsem scale)
fit_regsem_0$optim_fit
```

    [1] 0.5187112

``` r
# loss function value in regsem (on regsem scale)
fit_regsem_0$fit 
```

    [1] 0.5187112

The objective function value and the loss function value are the same because of `lambda = 0` (i.e., no penalty is implemented). We suspect that the loss function value in **regsem** is half of ML loss function value, which is justified by the following code

``` r
# chi-square value in lavaan
fitmeasures(fit_lavaan)["chisq"]
```

       chisq 
    312.2642 

``` r
# reproduce chi-square value using objective function value in regsem
2 * fit_regsem_0$optim_fit * nrow(lavaan::HolzingerSwineford1939)
```

    [1] 312.2642

``` r
# reproduce chi-square value using loss function value in regsem
2 * fit_regsem_0$fit * nrow(lavaan::HolzingerSwineford1939)
```

    [1] 312.2642

Now we try to understand the relationship between objective function value, loss function value, and lambda in **regsem**. The regularizer value in **regsem** can be obtained by

``` r
# regularizer value in regsem (on regsem scale)
regularizer_value <- fit_regsem$optim_fit - fit_regsem$fit
```

The absolute sum of penalized loadings is

``` r
sum_pen_coef <- 
  sum(abs(fit_regsem$coefficients[c("f -> x2", 
                                    "f -> x3", 
                                    "f -> x7", 
                                    "f -> x8", 
                                    "f -> x9")]))
```

Since the regularizer value should equal to the absolute sum of penalized loadings multiplied by lambda, then `regularizer_value / sum_pen_coef` should be `0.1`, the lambda we specified.

``` r
regularizer_value / sum_pen_coef
```

    [1] 0.09941301

By these examinations, we conclude that the scale of function value in **lslx** is **twice** of the scale of function value in **regsem**, so des the scale of lambda. Therefore, we can compare the objective function values made by **lslx** and **regsem**

``` r
# objective function value in lslx (on the same scale)
r6_lslx$extract_numerical_condition()["objective_value"]
```

    objective_value 
           1.246158 

``` r
# objective function value in regsem (on the same scale)
2 * fit_regsem$optim_fit
```

    [1] 1.249436

By the fact that the objective value made by **regsem** is larger than that made by **lslx**, it seems that the solution obtained by **regsem** can be further improved. The optimality of **lslx** solution can be checked by the sub-gradient of objective function

``` r
r6_lslx$extract_objective_gradient()
```

                       [,1]
    x1<-1|G    0.000000e+00
    x2<-1|G    0.000000e+00
    x3<-1|G    0.000000e+00
    x4<-1|G    0.000000e+00
    x5<-1|G    0.000000e+00
    x6<-1|G    0.000000e+00
    x7<-1|G    0.000000e+00
    x8<-1|G    0.000000e+00
    x9<-1|G    0.000000e+00
    x1<-f1|G   2.165399e-02
    x2<-f1|G   0.000000e+00
    x3<-f1|G   0.000000e+00
    x4<-f1|G  -2.379782e-04
    x5<-f1|G  -1.557863e-04
    x6<-f1|G  -5.170913e-05
    x7<-f1|G   0.000000e+00
    x8<-f1|G   0.000000e+00
    x9<-f1|G   8.651358e-04
    f1<->f1|G -9.957826e-05
    x1<->x1|G  4.544209e-05
    x2<->x2|G  1.744168e-15
    x3<->x3|G -1.365978e-16
    x4<->x4|G -7.327752e-05
    x5<->x5|G  1.006417e-04
    x6<->x6|G -3.726896e-05
    x7<->x7|G -1.585970e-16
    x8<->x8|G -2.125534e-16
    x9<->x9|G  1.554042e-04

Note that the largest element in the sub-gradient is for coefficient `x1<-f1|G`, which is a fixed coefficient and hence it is larger than our specified convergence criterion.

The coefficient estimates yielded by **lslx** and **regsem** are

``` r
r6_lslx$extract_coefficient()
```

      x1<-1|G   x2<-1|G   x3<-1|G   x4<-1|G   x5<-1|G   x6<-1|G   x7<-1|G   x8<-1|G   x9<-1|G  x1<-f1|G 
    4.9357697 6.0880399 2.2504153 3.0609081 4.3405316 2.1855719 4.1859021 5.5270764 5.3741233 1.0000000 
     x2<-f1|G  x3<-f1|G  x4<-f1|G  x5<-f1|G  x6<-f1|G  x7<-f1|G  x8<-f1|G  x9<-f1|G f1<->f1|G x1<->x1|G 
    0.0000000 0.0000000 1.9942417 2.2140809 1.8372766 0.0000000 0.0000000 0.1043953 0.2474084 1.1357065 
    x2<->x2|G x3<->x3|G x4<->x4|G x5<->x5|G x6<->x6|G x7<->x7|G x8<->x8|G x9<->x9|G 
    1.3818839 1.2749649 0.3666130 0.4468876 0.3612504 1.1832395 1.0220828 0.9919392 

``` r
fit_regsem$coefficients
```

      f -> x2 f -> x3 f -> x4 f -> x5 f -> x6 f -> x7 f -> x8 f -> x9 1 -> x1 1 -> x2 1 -> x3 1 -> x4
    1       0   0.004   1.903   2.106   1.755   0.006   0.009   0.146   4.936   6.088    2.25   3.062
      1 -> x5 1 -> x6 1 -> x7 1 -> x8 1 -> x9 x1 ~~ x1 x2 ~~ x2 x3 ~~ x3 x4 ~~ x4 x5 ~~ x5 x6 ~~ x6
    1    4.34   2.185   4.186   5.527   5.374    1.122    1.344    1.245    0.377    0.455    0.371
      x7 ~~ x7 x8 ~~ x8 x9 ~~ x9 f ~~ f
    1    1.171    1.016    0.969  0.266

We can observe that the estimates are slightly different.

As a double check, we also fit the model to data via **lsl**, a predecessor of **lslx**.

``` r
library(lsl)
lambda <- matrix(NA, 9, 1)
lambda[c(1, 4, 5, 6), 1] <- 1

rc_sem <- lslSEM()
rc_sem$input(raw = lavaan::HolzingerSwineford1939)
rc_sem$specify(pattern = list(lambda = lambda))
rc_sem$learn(penalty = list(type = "l1", 
                            gamma = 0.2, 
                            delta = 5), 
             variable = 7:15)
rc_sem$summarize(type = "overall")["dpl", "bic optimal"]
```

    [1] 1.246392

``` r
rc_sem$summarize(type = "individual")[, "bic optimal", drop = FALSE]
```

                bic optimal
    lambda[2,1]   0.0000000
    lambda[3,1]   0.0000000
    lambda[4,1]   1.9975814
    lambda[5,1]   2.2178113
    lambda[6,1]   1.8401894
    lambda[7,1]   0.0000000
    lambda[8,1]   0.0000000
    lambda[9,1]   0.1016986
    psi[1,1]      1.1394471
    psi[2,2]      1.3863898
    psi[3,3]      1.2791144
    psi[4,4]      0.3677369
    psi[5,5]      0.4481618
    psi[6,6]      0.3623881
    psi[7,7]      1.1870833
    psi[8,8]      1.0253894
    psi[9,9]      0.9955782
    phi[1,1]      0.2474556
    nu[1,1]       4.9357697
    nu[2,1]       6.0880399
    nu[3,1]       2.2504153
    nu[4,1]       3.0609081
    nu[5,1]       4.3405316
    nu[6,1]       2.1855719
    nu[7,1]       4.1859021
    nu[8,1]       5.5270764
    nu[9,1]       5.3741233

The objective function value and estimate made by **lsl** is consistent with **lslx**, which further supports that solution of **regsem** can be further improved. Note that **lslx** and **lsl** implement totally different type of algorithms for optimization. The obtained consistent result cannot be attributed to the similarity of algorithms.

### Computation Time

Finally, we compare the computation time of **lslx** and **regsem**

``` r
fun_lslx <- function() {
model_lslx <-
"f1 :=> fix(1) * x1 + pen() * x2 + pen() * x3 + x4 + x5 + x6 + pen() * x7 + pen() * x8 + pen() * x9"
r6_lslx <- lslx$new(model = model_lslx,
                    data = lavaan::HolzingerSwineford1939,
                    verbose = FALSE)
r6_lslx$fit_lasso(lambda_grid = 0.2, verbose = FALSE)
}

fun_regsem <- function() {
model_lavaan <- 
"f =~ 1 * x1 + l1 * x2 + l2 * x3 + l3 * x4 + l4 * x5 + l5 * x6 + l6 * x7 + l7 * x8 + l8 * x9"
fit_lavaan <- lavaan::cfa(model = model_lavaan, 
                          data = lavaan::HolzingerSwineford1939, 
                          meanstructure = TRUE)
fit_regsem <- regsem(fit_lavaan, lambda = 0.1, type = "lasso",
                     pars_pen = c("l1", "l2", "l6", "l7", "l8"),
                     solver = TRUE, quasi = TRUE,
                     solver.maxit = 100, line.search = TRUE)
}

microbenchmark::microbenchmark(
  fun_lslx(),
  fun_regsem(),
  times = 10)
```

    Unit: milliseconds
             expr      min       lq      mean   median       uq        max neval
       fun_lslx() 15.10807 15.45935  19.60642 16.25128 19.19613   44.40554    10
     fun_regsem() 77.90568 78.81560 216.14107 81.13811 83.00181 1433.47623    10

We can see that **lslx** is about 20 times faster than **regsem**.
