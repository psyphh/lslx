## ----comment = "", message = FALSE, setup, include=FALSE------------------------------------------
options(digits = 3)
options(width = 100)

## ----comment = "", message = FALSE----------------------------------------------------------------
model <-
'
visual  :=> x1 + x2 + x3
textual :=> x4 + x5 + x6
speed   :=> x7 + x8 + x9
visual  :~> x4 + x5 + x6 + x7 + x8 + x9 
textual :~> x1 + x2 + x3 + x7 + x8 + x9 
speed   :~> x1 + x2 + x3 + x4 + x5 + x6 
visual  <=> fix(1) * visual
textual <=> fix(1) * textual
speed   <=> fix(1) * speed
'

## ----comment = "", message = FALSE----------------------------------------------------------------
library(lslx)
r6_lslx <- lslx$new(model = model,
                    data = lavaan::HolzingerSwineford1939)

## ----comment = "", message = FALSE----------------------------------------------------------------
r6_lslx$extract_specification()

## ----comment = "", message = FALSE----------------------------------------------------------------
r6_lslx$fit(penalty_method = "mcp",
            lambda_grid = seq(.01, .30, .01),
            delta_grid = c(5, 10))

## ----comment = "", message = FALSE, fig.width = 24, fig.height = 14-------------------------------
r6_lslx$summarize(selector = "bic")

## ----comment = "", message = FALSE, fig.width = 8, fig.height = 4, dpi=200, out.width=600, out.height=300----
r6_lslx$plot_numerical_condition()

## ----comment = "", message = FALSE, fig.width = 8, fig.height = 4, dpi=200, out.width=600, out.height=300----
r6_lslx$plot_information_criterion()

## ----comment = "", message = FALSE, fig.width = 8, fig.height = 4, dpi=200, out.width=600, out.height=300----
r6_lslx$plot_fit_indice()

## ----comment = "", message = FALSE, fig.width = 8, fig.height = 4, dpi=200, out.width=600, out.height=300----
r6_lslx$plot_coefficient(block = "y<-f")

## ----comment = "", message = FALSE, fig.width = 8, fig.height = 4, dpi=300, out.width=600, out.height=300----
r6_lslx$extract_coefficient_matrice(selector = "bic", block = "y<-f")

## ----comment = "", message = FALSE, fig.width = 8, fig.height = 4, dpi=300, out.width=600, out.height=300----
r6_lslx$extract_implied_cov(selector = "bic")

## ----comment = "", message = FALSE, fig.width = 8, fig.height = 4, dpi=300, out.width=600, out.height=300----
r6_lslx$extract_residual_cov(selector = "bic")

