lslx$set("public",
         "initialize",
         function(model,
                  data,
                  group_variable,
                  reference_group,
                  weight_variable,
                  auxiliary_variable,
                  sample_cov,
                  sample_mean,
                  sample_size,
                  verbose = TRUE) {
           if (missing(model)) {
             stop("Argument 'model' cannot be empty.")
           }
           if (missing(data) & missing(sample_cov)) {
             stop("Argument 'data' and 'sample_cov' cannot be both empty.")
           } else if (!missing(data)) {
             if (!is.data.frame(data)) {
               if (!(is.matrix(data) & is.numeric(data))) {
                 stop("Argument 'data' must be a 'data.frame' or a 'matrix'.")
               } else {
                 data <- as.data.frame(data)
               }
             }
             if (missing(group_variable)) {
               name_group <- "G"
             } else {
               if (length(group_variable) > 1) {
                 stop("Argument `group_variable` can be only of length one.")
               }
               if (!(group_variable %in% colnames(data))) {
                 stop("Argument 'group_variable' is not recognized.")
               }
               name_group <-
                 sort(levels(factor(getElement(
                   data, group_variable
                 ))))
             }
             if (missing(weight_variable)) {
               
             } else {
               if (length(weight_variable) > 1) {
                 stop("Argument `weight_variable` can be only of length one.")
               }
               if (!(weight_variable %in% colnames(data))) {
                 stop("Argument 'weight_variable' is not recognized.")
               }
             }
             if (missing(auxiliary_variable)) {
               
             } else {
               if (!all(auxiliary_variable %in% colnames(data))) {
                 stop("Argument 'auxiliary_variable' is not recognized.")
               }
             }
           } else {
             if (!is.matrix(sample_cov) & !is.list(sample_cov)) {
               stop(
                 "Argument 'sample_cov' must be a 'matrix' (for single group analysis)",
                 " or a 'list' of 'matrix' (for multiple group analysis)."
               )
             }
             if (is.matrix(sample_cov)) {
               sample_cov <- list(sample_cov)
             }
             if (is.null(names(sample_cov))) {
               if (length(sample_cov) > 1) {
                 name_group <- paste0("G", 1:length(sample_cov))
               } else {
                 name_group <- "G"
               }
               names(sample_cov) <- name_group
               if (verbose) {
                 cat(
                   "NOTE: Because argument 'sample_cov' doesn't contain group name(s),",
                   "default group name(s) is created.\n"
                 )
               }
             } else {
               name_group <- names(sample_cov)
             }
             if (!missing(group_variable)) {
               stop("Argument 'group_variable' is unnecessary under moment initialization.")
             }
             if (!missing(weight_variable)) {
               stop("Argument 'weight_variable' is unnecessary under moment initialization.")
             }
             if (!missing(auxiliary_variable)) {
               stop("Argument 'auxiliary_variable' is unnecessary under moment initialization.")
             }
           }
           if (any(grepl(pattern = "/|\\||@",
                         x = name_group))) {
             stop(
               "Name(s) of group(s) cannot contain '/', '|', and '@'.",
               "\n  Please change the name(s) of group(s) in the specified data source."
             )
           }
           
           if (missing(reference_group)) {
             reference_group <- NA_character_
           } else {
             if (length(name_group) == 1L) {
               stop("Argument 'reference_group' is unnecessary for single group analysis.")
             } else {
               if (!(reference_group %in% name_group)) {
                 stop(
                   "Argument 'reference_group' is not recognized.",
                   "\n  Group name(s) currently recognized by 'lslx' is ",
                   do.call(paste, as.list(name_group)),
                   " (possibly automatically created).",
                   "\n  Specified 'reference_group' is ",
                   reference_group,
                   "."
                 )
               }
             }
           }
           private$model <-
             lslxModel$new(model = model,
                           name_group = name_group,
                           reference_group = reference_group)
           if (!missing(data)) {
             if (!all(private$model$name_response %in% colnames(data))) {
               stop(
                 "Some response variable in 'model' cannot be found in 'data'.",
                 "\n  Response variables specified by 'model' are ",
                 do.call(paste, as.list(private$model$name_response)),
                 ".",
                 "\n  Column names of 'data' are ",
                 do.call(paste, as.list(colnames(data))),
                 "."
               )
             } else {
               if (missing(group_variable)) {
                 response <-
                   list(data[, private$model$name_response, drop = FALSE])
                 names(response) <- name_group
                 if (missing(weight_variable)) {
                   weight <- list(data.frame(weight = rep(1, nrow(data))))
                 } else {
                   weight <- list(data[, weight_variable, drop = FALSE])
                 }
                 names(weight) <- name_group
                 if (missing(auxiliary_variable)) {
                   auxiliary <- list()
                 } else {
                   auxiliary_variable <-
                     setdiff(x = auxiliary_variable, 
                             y = private$model$name_response)
                   if (length(auxiliary_variable) > 0) {
                     auxiliary <- list(data[, auxiliary_variable, drop = FALSE])
                   } else {
                     auxiliary <- list()
                   }
                 }
               } else {
                 data <-
                   data[order(as.character(getElement(data, group_variable))),]
                 data[, group_variable] <-
                   as.character(getElement(data, group_variable))
                 response <-
                   split(data[, private$model$name_response, drop = FALSE],
                         getElement(data, group_variable))
                 if (missing(weight_variable)) {
                   weight <- split(data.frame(weight = rep(1, nrow(data))),
                                   getElement(data, group_variable))
                 } else {
                   weight <-
                     split(data[, weight_variable, drop = FALSE],
                           getElement(data, group_variable))
                 }
                 if (missing(auxiliary_variable)) {
                   auxiliary <- list()
                 } else {
                   auxiliary <-
                     split(data[, auxiliary_variable, drop = FALSE],
                           getElement(data, group_variable))
                 }
               }
             }
             if (!all(sapply(X = response, 
                             FUN = function(response_i) {
                               sapply(X = response_i,
                                      FUN = function(response_ij) {
                                        return(is.numeric(response_ij))
                                      })
                             }))) {
               stop("Response variable(s) cannot contain non-numeric variables.")
             }
             if (length(auxiliary) > 0) {
               if (!all(sapply(X = auxiliary, 
                               FUN = function(auxiliary_i) {
                                 sapply(X = auxiliary_i,
                                        FUN = function(auxiliary_ij) {
                                          return(is.numeric(auxiliary_ij))
                                        })
                               }))) {
                 stop("Auxiliary variable(s) cannot contain non-numeric variables.")
               }
             }
             if (any(sapply(
               X = weight,
               FUN = function(weight_i) {
                 return(any(weight_i < 0))
               }
             ))) {
               stop("Weight variable cannot contain negative value.")
             }
             private$data <- lslxData$new(response = response,
                                          weight = weight,
                                          auxiliary = auxiliary)
             if (verbose) {
               cat("An 'lslx' R6 class is initialized via 'data'.\n")
             }
           } else {
             if (!all(private$model$name_response %in% colnames(sample_cov[[1]]))) {
               stop(
                 "Some response variable in 'model' cannot be found in 'sample_cov'.",
                 "\n  Response variables specified by 'model' are ",
                 do.call(paste, as.list(private$model$name_response)),
                 ".",
                 "\n  Column names of 'sample_cov' are ",
                 do.call(paste, as.list(colnames(sample_cov[[1]]))),
                 ".",
                 "\n  Row names of 'sample_cov' are ",
                 do.call(paste, as.list(rownames(sample_cov[[1]]))),
                 "."
               )
             } else {
               sample_cov <-
                 lapply(
                   X = sample_cov,
                   FUN = function(sample_cov_i) {
                     sample_cov_i <-
                       sample_cov_i[private$model$name_response,
                                    private$model$name_response,
                                    drop = FALSE]
                     return(sample_cov_i)
                   }
                 )
               names(sample_cov) <- name_group
             }
             if (missing(sample_mean)) {
               sample_mean <-
                 lapply(
                   X = name_group,
                   FUN = function(i) {
                     sample_mean_i <- rep(0, ncol(sample_cov[[1]]))
                     names(sample_mean_i) <-
                       colnames(sample_cov[[1]])
                     return(sample_mean_i)
                   }
                 )
               names(sample_mean) <- name_group
               if (verbose) {
                 cat(
                   "NOTE: Because argument 'sample_mean' is missing,",
                   "default 'sample_mean' is created.\n"
                 )
               }
             } else {
               if (!is.numeric(sample_mean) & !is.list(sample_mean)) {
                 stop("Argument 'sample_mean' must be 'numeric' or 'list' of 'numeric'.")
               }
               if (is.numeric(sample_mean)) {
                 sample_mean <- list(sample_mean)
               }
               if (!all(private$model$name_response %in% names(sample_mean[[1]]))) {
                 stop(
                   "Some response variable in 'model' cannot be found in 'sample_mean'.",
                   "\n  Response variables specified by 'model' are ",
                   do.call(paste, as.list(private$model$name_response)),
                   ".",
                   "\n  Column names of 'sample_mean' are ",
                   do.call(paste, as.list(colnames(
                     sample_cov[[1]]
                   ))),
                   "."
                 )
               } else {
                 sample_mean <-
                   lapply(
                     X = sample_mean,
                     FUN = function(sample_mean_i) {
                       sample_mean_i <-
                         sample_mean_i[private$model$name_response]
                       return(sample_mean_i)
                     }
                   )
               }
               if (length(sample_mean) != length(sample_cov)) {
                 stop(
                   "The length of argument 'sample_mean' doesn't match the length of 'sample_cov'.",
                   "\n  The length of 'sample_mean' is ",
                   length(sample_mean),
                   ".",
                   "\n  The length of 'sample_cov' is ",
                   length(sample_cov),
                   "."
                 )
               }
               if (is.null(names(sample_mean))) {
                 names(sample_mean) <- name_group
               } else {
                 if (!all(name_group %in% names(sample_mean))) {
                   stop(
                     "Argument 'sample_mean' contains unrecognized group name(s).",
                     "\n  Group name(s) currently recognized by 'lslx' is ",
                     do.call(paste, as.list(private$model$name_group)),
                     " (possibly automatically created).",
                     "\n  Group name(s) specified in 'sample_mean' is ",
                     do.call(paste, as.list(names(
                       sample_mean
                     ))),
                     "."
                   )
                 } else{
                   sample_mean <- sample_mean[name_group]
                 }
               }
             }
             if (missing(sample_size)) {
               stop("Argument 'sample_size' cannot be empty if 'sample_cov' is used.")
             } else {
               if (!is.numeric(sample_size) & !is.list(sample_size)) {
                 stop("Argument 'sample_size' must be a 'numeric' or a 'list' of 'numeric'.")
               }
               if (is.numeric(sample_size)) {
                 sample_size <- as.list(sample_size)
               }
               if (length(sample_size) != length(sample_cov)) {
                 stop(
                   "The length of argument 'sample_size' doesn't match the length of 'sample_cov'.",
                   "\n  The length of 'sample_size' is ",
                   length(sample_size),
                   ".",
                   "\n  The length of 'sample_cov' is ",
                   length(sample_cov),
                   "."
                 )
               }
               if (is.null(names(sample_size))) {
                 names(sample_size) <- name_group
               } else {
                 if (!all(name_group %in% names(sample_size))) {
                   stop(
                     "Argument 'sample_size' contains unrecognized group name(s).",
                     "\n  Group name(s) currently recognized by 'lslx' is ",
                     do.call(paste, as.list(private$model$name_group)),
                     " (possibly automatically created).",
                     "\n  Group name(s) specified in 'sample_size' is ",
                     do.call(paste, as.list(names(
                       sample_size
                     ))),
                     "."
                   )
                 } else{
                   sample_size <- sample_size[name_group]
                 }
               }
             }
             private$data <-
               lslxData$new(
                 sample_cov = sample_cov,
                 sample_mean = sample_mean,
                 sample_size = sample_size
               )
             if (verbose) {
               cat("An 'lslx' R6 class is initialized via 'sample_cov'. \n")
             }
           }
           private$fitting <- NULL
           if (verbose) {
             cat("  Response Variable(s):",
                 private$model$name_response,
                 "\n")
             if (length(private$model$name_factor) > 0) {
               cat("  Latent Factor(s):",
                   private$model$name_factor,
                   "\n") 
             }
             if (length(private$data$auxiliary) > 0) {
               cat("  Auxiliary Variable(s):",
                   colnames(private$data$auxiliary[[1]]),
                   "\n")               
             }
             if (length(private$model$name_group) > 1) {
               cat("  Group(s):",
                   private$model$name_group,
                   "\n")
               if (!is.na(private$model$reference_group)) {
                 cat("  Reference Group:",
                     private$model$reference_group,
                     "\n") 
               }
             }
             if (!is.na(private$model$reference_group)) {
               cat(
                 "NOTE:",
                 "Because",
                 private$model$reference_group,
                 "is set as reference,",
                 "coefficients in other groups actually represent increments from the reference.\n"
               )
             }
           }
         })
