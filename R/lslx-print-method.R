lslx$set("public",
         "print",
         function() {
           cat("To fit model to data, please use\n")
           cat("$fit(), fit_lasso(), fit_mcp(), or fit_none()")
           invisible(self)
         })
