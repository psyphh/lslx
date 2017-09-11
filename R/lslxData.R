lslxData <-
  R6::R6Class(
    classname = "lslxData",
    public = list(
      sample_data = "list",
      sample_cov = "list",
      sample_mean = "list",
      sample_size = "list"
    )
  )



lslxData$set("public",
             "initialize",
             function(sample_data,
                      sample_cov,
                      sample_mean,
                      sample_size) {
               if (!missing(sample_data)) {
                 self$sample_data <- sample_data
                 self$sample_size <-
                   lapply(X = sample_data, nrow)
               } else {
                 self$sample_data <- list()
                 self$sample_cov <- sample_cov
                 self$sample_mean <- sample_mean
                 self$sample_size <- sample_size
               }
             })
