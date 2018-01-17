## ----comment = "", message = FALSE, setup, include=FALSE------------------------------------------
options(digits = 3)
options(width = 100)

## ----comment = "", message = FALSE----------------------------------------------------------------
data <- lavaan::HolzingerSwineford1939
data$x5 <- ifelse(data$x1 <= quantile(data$x1, .3), NA, data$x5)
data$age <- data$ageyr + data$agemo/12
data$x9 <- ifelse(data$age <= quantile(data$age, .3), NA, data$x9)

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
                    data = data,
                    auxiliary_variable = c("ageyr", "agemo"))

## ----comment = "", message = FALSE----------------------------------------------------------------
r6_lslx$fit(penalty_method = "mcp",
            lambda_grid = seq(.01, .30, .01),
            delta_grid = c(5, 10))

## ----comment = "", message = FALSE, fig.width = 24, fig.height = 14-------------------------------
r6_lslx$summarize(selector = "bic")

