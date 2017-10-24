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
                 idc_use <-
                   lapply(X = self$pattern,
                          FUN = function(pattern_i) {
                            return(rowSums(pattern_i) < nrow(pattern_i))
                          })
                 self$response <- 
                   mapply(FUN = function(response_i,
                                         idc_use_i) {
                            return(response_i[idc_use_i, , drop = FALSE])
                          },
                          self$response,
                          idc_use,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE)
                 
                 self$weight <- 
                   mapply(FUN = function(weight_i,
                                         idc_use_i) {
                     weight_i <- weight_i[idc_use_i]
                     return(weight_i / sum(weight_i))
                   },
                   self$weight,
                   idc_use,
                   SIMPLIFY = FALSE,
                   USE.NAMES = TRUE)
                 
                 self$pattern <- 
                   mapply(FUN = function(pattern_i,
                                         idc_use_i) {
                     return(pattern_i[idc_use_i, ])
                   },
                   self$pattern,
                   idc_use,
                   SIMPLIFY = FALSE,
                   USE.NAMES = TRUE)
                 
                 if (length(auxiliary) > 0) {
                   self$auxiliary <- 
                     mapply(FUN = function(auxiliary_i,
                                           idc_use_i) {
                       return(auxiliary_i[idc_use_i, ])
                     },
                     self$auxiliary,
                     idc_use,
                     SIMPLIFY = FALSE,
                     USE.NAMES = TRUE)
                 } 
                 
                 self$sample_cov <- list()
                 self$sample_mean <- list()
                 self$sample_size <-
                   lapply(X = response, 
                          FUN = nrow)
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
