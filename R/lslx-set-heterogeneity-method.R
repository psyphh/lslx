lslx$set("private",
         "set_heterogeneity",
         function(block,
                  group,
                  type,
                  verbose = TRUE) {
           
           if (is.na(private$model$reference_group)){
             stop("Reference group is not specified.",
                  "\n The 'set heterogeneity method' can be used only if reference group is specified")
           }
           if (group==private$model$reference_group){
             stop("Argument 'group' is set as the reference group, please select another group.",
                  "\n  Group name(s) currently recognized by 'lslx' is ",
                  do.call(paste, as.list(private$model$name_group)),
                  ".",
                  "\n  Group name specified in 'group' is ",
                  group,
                  ".")
           }else if (!(group %in% private$model$name_group)) {
             stop(
               "Argument 'group' contains unknown group name.",
               "\n  Group name(s) currently recognized by 'lslx' is ",
               do.call(paste, as.list(private$model$name_group)),
               ".",
               "\n  Group name specified in 'group' is ",
               group,
               "."
             )
           }
           relation<-private$model$specification[(private$model$specification$block==block)&
                                                   (private$model$specification$type!="fixed")&
                                                   private$model$specification$reference,"relation"]
            if (length(relation)==0) {
             stop("No valid relation ",
                  block,
                  " under group ",
                  group,
                  " is found. Please check the settings.")
            } else {
              name<-paste0(relation,"|",group)
              private$set_coefficient(
                name = name,
                type = type,
                verbose = verbose
              )
            }
           
         }
)


lslx$set("public",
         "free_heterogeneity",
         function(block,
                  group,
                  verbose = TRUE) {
           private$set_heterogeneity(block = block,
                                     group = group,
                                     type = "free",
                                     verbose = verbose)
         }
)


lslx$set("public",
         "fix_heterogeneity",
         function(block,
                  group,
                  verbose = TRUE) {
           private$set_heterogeneity(block = block,
                                     group = group,
                                     type = "fixed",
                                     verbose = verbose)
         }
)


lslx$set("public",
         "penalize_heterogeneity",
         function(block,
                  group,
                  verbose = TRUE) {
           private$set_heterogeneity(block = block,
                                     group = group,
                                     type = "pen",
                                     verbose = verbose)
         }
)
