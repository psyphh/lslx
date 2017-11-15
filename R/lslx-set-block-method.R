lslx$set("private",
         "set_block",
         function(block,
                  group,
                  type,
                  verbose = TRUE) {
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
           }
           relation<-private$model$specification[(private$model$specification$block %in% block)&
                                                   (private$model$specification$group %in% group),
                                                 "relation"]
           
           if (length(relation)==0) {
             stop("No valid relation ",
                  block,
                  " under group ",
                  group,
                  " is found. Please check the settings.")
           } else {
             name<-unique(paste0(expand.grid(relation,group)[,1],"|",expand.grid(relation,group)[,2]))
             private$set_coefficient(
               name = name,
               type = type,
               verbose = verbose
             )
           }
         }
)


lslx$set("public",
         "free_block",
         function(block,
                  group = "default",
                  verbose = TRUE) {
           private$set_block(block = block,
                             group = group,
                             type = "free",
                             verbose = verbose)
         }
)


lslx$set("public",
         "fix_block",
         function(block,
                  group = "default",
                  verbose = TRUE) {
           private$set_block(block = block,
                             group = group,
                             type = "fixed",
                             verbose = verbose)
         }
)


lslx$set("public",
         "penalize_block",
         function(block,
                  group = "default",
                  verbose = TRUE) {
           private$set_block(block = block,
                             group = group,
                             type = "pen",
                             verbose = verbose)
         }
)
