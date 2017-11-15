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
           
           if (group == "default") {
             group <-  private$model$name_group
           } else if (!(group %in% private$model$name_group)) {
             stop(
               "Argument 'group' contains unknown group name.",
               "\n  Group name(s) currently recognized by 'lslx' is ",
               do.call(paste, as.list(private$model$name_group)),
               ".",
               "\n  Group name specified in 'group' is ",
               group,
               "."
             )
           } else {}
           name <- paste0(expand.grid(left,"<-",right)[,1],
                          expand.grid(left,"<-",right)[,2],
                          expand.grid(left,"<-",right)[,3],
                          "|",group)
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
                  group = "default",
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
                  group = "default",
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
                  group = "default",
                  verbose = TRUE) {
           private$set_directed(
             left = left,
             right = right,
             group = group,
             type = "pen",
             verbose = verbose)
         }
)
