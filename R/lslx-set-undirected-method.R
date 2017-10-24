lslx$set("private",
         "set_undirected",
         function(both,
                  group,
                  type,
                  verbose = TRUE) {
           if (is.na(private$model$reference_group)) {
             name <- combn(both,2,function(x) paste0(x[1],"<->",x[2]))
           } else if (missing(group)) {
             stop("Argument 'group' must be given")
           } else { 
             name <- paste0(combn(both,2,function(x) paste0(x[1],"<->",x[2])),"|",group)
           }
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
                  group,
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
                  group,
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
                  group,
                  verbose = TRUE) {
           private$set_undirected(both = both,
                                  group = group,
                                  type = "pen",
                                  verbose = verbose)
         }
)
