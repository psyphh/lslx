## ----comment = "", message = FALSE, setup, include=FALSE------------------------------------------
options(digits = 3)
options(width = 100)

## ----comment = "", message = FALSE----------------------------------------------------------------
model <-
'
fix(1) * x1 + x2 + x3      <=: ind60
fix(1) * y1 + y2 + y3 + y4 <=: dem60
fix(1) * y5 + y6 + y7 + y8 <=: dem65
dem60 <= ind60
dem65 <= ind60 + dem60
'

## ----comment = "", message = FALSE----------------------------------------------------------------
library(lslx)
r6_lslx <- lslx$new(model = model,
                    sample_cov = cov(lavaan::PoliticalDemocracy),
                    sample_size = nrow(lavaan::PoliticalDemocracy))

## ----comment = "", message = FALSE----------------------------------------------------------------
r6_lslx$penalize_coefficient(name = c("y1<->y5",
                                      "y2<->y4",
                                      "y2<->y6",
                                      "y3<->y7",
                                      "y4<->y8",
                                      "y6<->y8"))

## ----comment = "", message = FALSE----------------------------------------------------------------
r6_lslx$fit_mcp(lambda_grid = seq(.01, .30, .01),
                 delta_grid = Inf)

## ----comment = "", message = FALSE, fig.width = 24, fig.height = 14-------------------------------
r6_lslx$summarize(selector = "aic")

## ----comment = "", message = FALSE, fig.width = 8, fig.height = 4, dpi=300, out.width=600, out.height=300----
r6_lslx$extract_coefficient(selector = "bic")

## ----comment = "", message = FALSE, fig.width = 8, fig.height = 4, dpi=300, out.width=600, out.height=300----
diag(r6_lslx$extract_coefficient_acov(selector = "bic"))

