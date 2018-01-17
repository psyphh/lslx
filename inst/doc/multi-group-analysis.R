## ----comment = "", message = FALSE, setup, include=FALSE------------------------------------------
options(digits = 3)
options(width = 100)

## ----comment = "", message = FALSE----------------------------------------------------------------
model <-
'
visual  :=> fix(1) * x1 + x2 + x3 
textual :=> fix(1) * x4 + x5 + x6 
speed   :=> fix(1) * x7 + x8 + x9 
'

## ----comment = "", message = FALSE----------------------------------------------------------------
library(lslx)
r6_lslx <- lslx$new(model = model,
                     data = lavaan::HolzingerSwineford1939,
                     group_variable = "school",
                     reference_group = "Pasteur")

## ----comment = "", message = FALSE----------------------------------------------------------------
r6_lslx$penalize_heterogeneity(block = "y<-f", group = "Grant-White")
r6_lslx$penalize_heterogeneity(block = "y<-1", group = "Grant-White")

## ----comment = "", message = FALSE----------------------------------------------------------------
r6_lslx$free_directed(left = c("visual", "textual", "speed"),
                      right = "1",
                      group = "Grant-White")

## ----comment = "", message = FALSE----------------------------------------------------------------
r6_lslx$fit_lasso(lambda_grid = seq(.01, .30, .01))

## ----comment = "", message = FALSE, fig.width = 24, fig.height = 14-------------------------------
r6_lslx$summarize(selector = "bic")

