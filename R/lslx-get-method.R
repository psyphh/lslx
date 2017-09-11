lslx$set("public",
         "get_model",
         function() {
           return(private$model$clone(deep = TRUE))
         })

lslx$set("public",
         "get_data",
         function() {
           return(private$data$clone(deep = TRUE))
         })

lslx$set("public",
         "get_fitting",
         function() {
           return(private$fitting$clone(deep = TRUE))
         })
