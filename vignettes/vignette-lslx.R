### R code from vignette source 'vignette-lslx.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: preliminaries
###################################################
library(lslx)
library(ggplot2)
library(dplyr)
library(knitr)
library(tidyr)
options(prompt = "R> ", continue = "+  ", width = 80, useFancyQuotes = FALSE)
options(digits = 3L)


###################################################
### code chunk number 2: vignette-lslx.Rnw:352-354
###################################################
model_reg <- "y <= x1 + x2
              y <~ x3 + x4"


###################################################
### code chunk number 3: vignette-lslx.Rnw:367-369
###################################################
model_reg <- "x1 + x2 => y
              x3 + x4 ~> y"


###################################################
### code chunk number 4: vignette-lslx.Rnw:372-373
###################################################
model_reg <- "y <= x1 + x2 + pen() * x3 + pen() * x4"


###################################################
### code chunk number 5: vignette-lslx.Rnw:376-377
###################################################
model_reg <- "y <~ free() * x1 + free() * x2 + x3 + x4"


###################################################
### code chunk number 6: vignette-lslx.Rnw:467-476
###################################################
model_fa <- "visual  :=> x1 + x2 + x3
             textual :=> x4 + x5 + x6
             speed   :=> x7 + x8 + x9
             visual  :~> x4 + x5 + x6 + x7 + x8 + x9
             textual :~> x1 + x2 + x3 + x7 + x8 + x9
             speed   :~> x1 + x2 + x3 + x4 + x5 + x6
             visual  <=> fix(1) * visual
             textual <=> fix(1) * textual
             speed   <=> fix(1) * speed"


###################################################
### code chunk number 7: vignette-lslx.Rnw:481-490
###################################################
model_fa_lavaan <- "visual  =~ x1 + x2 + x3
                    textual =~ x4 + x5 + x6
                    speed   =~ x7 + x8 + x9
                    pen() * visual  =~ x4 + x5 + x6 + x7 + x8 + x9
                    pen() * textual =~ x1 + x2 + x3 + x7 + x8 + x9
                    pen() * speed   =~ x1 + x2 + x3 + x4 + x5 + x6
                    visual  ~~ 1 * visual
                    textual ~~ 1 * textual
                    speed   ~~ 1 * speed"


###################################################
### code chunk number 8: vignette-lslx.Rnw:497-499
###################################################
lslx_fa <- lslx$new(model = model_fa, 
  data = lavaan::HolzingerSwineford1939)


###################################################
### code chunk number 9: vignette-lslx.Rnw:504-506
###################################################
lslx_fa$fit(penalty_method = "mcp",
  lambda_grid = seq(.01, .60, .01), delta_grid = c(1.5, 3.0, Inf))


###################################################
### code chunk number 10: vignette-lslx.Rnw:511-512
###################################################
lslx_fa$summarize(selector = "bic", interval = FALSE)


###################################################
### code chunk number 11: vignette-lslx.Rnw:603-604
###################################################
lslx_fa$plot_numerical_condition()


###################################################
### code chunk number 12: vignette-lslx.Rnw:613-614 (eval = FALSE)
###################################################
## lslx_fa$plot_numerical_condition()


###################################################
### code chunk number 13: vignette-lslx.Rnw:618-619 (eval = FALSE)
###################################################
## lslx_fa$plot_coefficient(block = "y<-f")


###################################################
### code chunk number 14: vignette-lslx.Rnw:626-627
###################################################
lslx_fa$plot_coefficient(block = "y<-f")


###################################################
### code chunk number 15: vignette-lslx.Rnw:645-648
###################################################
moment_jacobian <- lslx_fa$extract_moment_jacobian(
  selector = "bic", type = "effective")
min(svd(moment_jacobian)$d)


###################################################
### code chunk number 16: vignette-lslx.Rnw:685-693
###################################################
data_miss <- lavaan::HolzingerSwineford1939
data_miss$x5 <- ifelse(
  test = data_miss$x1 <= quantile(data_miss$x1, .3), 
  yes = NA, no = data_miss$x5)
data_miss$age <- data_miss$ageyr + data_miss$agemo / 12
data_miss$x9 <- ifelse(
  test = data_miss$age <= quantile(data_miss$age, .3), 
  yes = NA, no = data_miss$x9)


###################################################
### code chunk number 17: vignette-lslx.Rnw:696-704
###################################################
model_miss <- "x1 + x2 + x3 <=: visual
               x4 + x5 + x6 <=: textual
               x7 + x8 + x9 <=: speed
               visual  <=> 1 * visual
               textual <=> 1 * textual
               speed   <=> 1 * speed"
lslx_miss <- lslx$new(model = model_miss, data = data_miss,
  auxiliary_variable = c("ageyr", "agemo"), verbose = FALSE)


###################################################
### code chunk number 18: vignette-lslx.Rnw:707-708
###################################################
lslx_miss$penalize_block(block = "y<->y", type = "fixed", verbose = FALSE)


###################################################
### code chunk number 19: vignette-lslx.Rnw:711-712
###################################################
lslx_miss$fit_lasso(verbose = FALSE)


###################################################
### code chunk number 20: vignette-lslx.Rnw:715-716
###################################################
lslx_miss$summarize(selector = "raic", style = "minimal")


###################################################
### code chunk number 21: vignette-lslx.Rnw:719-720
###################################################
lslx_miss$extract_coefficient_matrix(selector = "raic", block = "y<->y")


###################################################
### code chunk number 22: vignette-lslx.Rnw:723-724
###################################################
lslx_miss$extract_fit_index(selector = "raic")


###################################################
### code chunk number 23: vignette-lslx.Rnw:762-768
###################################################
model_mgfa <- "1 * x1 + x2 + x3 <=: visual 
               1 * x4 + x5 + x6 <=: textual
               1 * x7 + x8 + x9 <=: speed"
lslx_mgfa <- lslx$new(model = model_mgfa,
  data = lavaan::HolzingerSwineford1939, group_variable = "school",
  reference_group = "Pasteur", verbose = FALSE)


###################################################
### code chunk number 24: vignette-lslx.Rnw:773-776
###################################################
model_mgfa <- "c(fix(0), fix(1)) * x1 + x2 + x3 <=: visual 
               c(fix(0), fix(1)) * x4 + x5 + x6 <=: textual
               c(fix(0), fix(1)) * x7 + x8 + x9 <=: speed"


###################################################
### code chunk number 25: vignette-lslx.Rnw:781-783
###################################################
lslx_mgfa$penalize_heterogeneity(block = c("y<-f", "y<-1"), 
  group = "Grant-White", verbose = FALSE)


###################################################
### code chunk number 26: vignette-lslx.Rnw:786-788
###################################################
lslx_mgfa$free_block(block = "f<-1", 
  group = "Grant-White", verbose = FALSE)


###################################################
### code chunk number 27: vignette-lslx.Rnw:791-792
###################################################
lslx_mgfa$fit_mcp(verbose = FALSE)


###################################################
### code chunk number 28: vignette-lslx.Rnw:795-801
###################################################
loading <- lslx_mgfa$extract_coefficient_matrix(
  selector = "hbic", block = "y<-f")
intercept <- lslx_mgfa$extract_coefficient_matrix(
  selector = "hbic", block = "y<-1")
loading$"Grant-White" - loading$"Pasteur"
t(intercept$"Grant-White" - intercept$"Pasteur")


