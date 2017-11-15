lslx$set("private",
         "set_undirected",
         function(both,
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
           } else {}
           name <- paste0(
             expand.grid(combn(both,2,function(x) paste0(x[1],"<->",x[2])),group)[,1],
             "|",
             expand.grid(combn(both,2,function(x) paste0(x[1],"<->",x[2])),group)[,2]
           )
           private$set_coefficient(
             name = name,
             type = type,
             verbose = verbose
           )
         }
)


lslx$set("public",
         "free_undirected",
         function(both,
                  group = "default",
                  verbose = TRUE) {
           private$set_undirected(both = both,
                                  group = group,
                                  type = "free",
                                  verbose = verbose)
         }
)


lslx$set("public",
         "fix_undirected",
         function(both,
                  group = "default",
                  verbose = TRUE) {
           private$set_undirected(both = both,
                                  group = group,
                                  type = "fixed",
                                  verbose = verbose)
         }
)


lslx$set("public",
         "penalize_undirected",
         function(both,
                  group = "default",
                  verbose = TRUE) {
           private$set_undirected(both = both,
                                  group = group,
                                  type = "pen",
                                  verbose = verbose)
         }
)
