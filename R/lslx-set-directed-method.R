lslx$set("private",
         "set_directed",
         function(left,
                  right,
                  group,
                  type,
                  verbose = TRUE) {
           if (missing(left)) {
             stop("Argument 'left' must be given.")
           } else if(missing(right)) {
             stop("Argument 'right' must be given.")
           } else {}
           
           combination <- expand.grid(left,"<-",right)
           combination <- apply(combination,c(1,2), as.character)
           
           if (is.na(private$model$reference_group)) {
             name <- paste0(combination[,1],combination[,2],combination[,3])
           } else if (missing(group)) {
             stop("Argument 'group' must be given")
           } else { 
             name <- paste0(combination[,1],combination[,2],combination[,3],"|",group)
           }
           
           private$set_coefficient(
             name = name,
             type = type,
             verbose = verbose
           )
         }
)


lslx$set("public",
         "free_directed",
         function(left,
                  right,
                  group,
                  verbose = TRUE) {
           private$set_directed(
             left = left,
             right = right,
             group = group,
             type = "free",
             verbose = verbose)
         }
)


lslx$set("public",
         "fix_directed",
         function(left,
                  right,
                  group,
                  verbose = TRUE) {
           private$set_directed(
             left = left,
             right = right,
             group = group,
             type = "fixed",
             verbose = verbose)
         }
)


lslx$set("public",
         "penalize_directed",
         function(left,
                  right,
                  group,
                  verbose = TRUE) {
           private$set_directed(
             left = left,
             right = right,
             group = group,
             type = "pen",
             verbose = verbose)
         }
)
