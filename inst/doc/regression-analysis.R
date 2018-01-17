## ----comment = "", message = FALSE, setup, include=FALSE------------------------------------------
options(digits = 3)
options(width = 100)

## ----comment = "", message = FALSE----------------------------------------------------------------
set.seed(9487)
x <- matrix(rnorm(2000), 200, 10)
colnames(x) <- paste0("x", 1:10)
y <- matrix(rnorm(200), 200, 1)
data <- data.frame(y, x)

## ----comment = "", message = FALSE----------------------------------------------------------------
model <-
'
y <= x1 + x2 + x3 + x4
y <~ x5 + x6 + x7 + x8 + x9 + x10
'

## ----comment = "", message = FALSE----------------------------------------------------------------
library(lslx)
r6_lslx <- lslx$new(model = model, data = data)

## ----comment = "", message = FALSE----------------------------------------------------------------
r6_lslx$fit(penalty_method = "lasso",
            lambda_grid = seq(.00, .30, .01))

## ----comment = "", message = FALSE, fig.width = 24, fig.height = 14-------------------------------
r6_lslx$summarize(selector = "aic")

## ----comment = "", message = FALSE, fig.width = 8, fig.height = 4, dpi=200, out.width=600, out.height=300----
r6_lslx$plot_numerical_condition()

## ----comment = "", message = FALSE, fig.width = 8, fig.height = 4, dpi=200, out.width=600, out.height=300----
r6_lslx$plot_information_criterion()

## ----comment = "", message = FALSE, fig.width = 8, fig.height = 4, dpi=200, out.width=600, out.height=300----
r6_lslx$plot_fit_indice()

## ----comment = "", message = FALSE, fig.width = 8, fig.height = 4, dpi=200, out.width=600, out.height=300----
r6_lslx$plot_coefficient(block = "y<-y")

## ----comment = "", message = FALSE, fig.width = 8, fig.height = 4, dpi=300, out.width=600, out.height=300----
r6_lslx$extract_coefficient(selector = "aic")

## ----comment = "", message = FALSE, fig.width = 8, fig.height = 4, dpi=300, out.width=600, out.height=300----
r6_lslx$extract_objective_gradient(selector = "aic")

