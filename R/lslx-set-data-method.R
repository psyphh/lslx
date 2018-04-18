## \code{$set_data()} reset the data field. ##
lslx$set("public",
         "set_data",
         function(data,
                  sample_cov,
                  sample_mean,
                  sample_size) {
           private$data <- 
             lslxData$new(data = data,
                          sample_cov = sample_cov,
                          sample_mean = sample_mean,
                          sample_size = sample_size,
                          group_variable = private$model$group_variable,
                          weight_variable = private$model$weight_variable,
                          auxiliary_variable = private$model$auxiliary_variable,
                          name_response = private$model$name_response,
                          name_group = private$model$name_group)
           private$fitting <- NULL
         })
