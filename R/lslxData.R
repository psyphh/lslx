lslxData <-
  R6::R6Class(
    classname = "lslxData",
    public = list(
      response = "list",
      pattern = "list",
      weight = "list",
      auxiliary = "list",
      sample_cov = "list",
      sample_mean = "list",
      sample_size = "list"
    )
  )



lslxData$set("public",
             "initialize",
             function(response,
                      weight,
                      auxiliary,
                      sample_cov,
                      sample_mean,
                      sample_size) {
               if (!missing(response)) {
                 self$response <- response
                 self$weight <- weight
                 self$auxiliary <- auxiliary
                 self$pattern <- 
                   lapply(X = self$response,
                          FUN = function(response_i) {
                            return(!is.na(response_i))
                          })
                 self$sample_cov <- list()
                 self$sample_mean <- list()
                 self$sample_size <- list()
               } else {
                 self$response <- list()
                 self$weight <- list()
                 self$auxiliary <- list()
                 self$pattern <- list()
                 self$sample_cov <- sample_cov
                 self$sample_mean <- sample_mean
                 self$sample_size <- sample_size
               }
             })
