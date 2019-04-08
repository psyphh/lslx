## define R6 class \code{lslxModel} to store model information. ##
lslxModel <-
  R6::R6Class(
    classname = "lslxModel",
    public = list(
      model = "character",
      numeric_variable = "character",
      ordered_variable = "character",
      weight_variable = "character",
      auxiliary_variable = "character",
      group_variable = "character",
      reference_group = "character",
      level_group = "character",
      name_response = "character",
      name_factor = "character",
      name_eta = "character",
      name_exogenous = "character",
      name_endogenous = "character",
      name_covariate = "character",
      nlevel_ordered = "numeric",
      specification = "data.frame"
    )
  )

## \code{$new()} initializes a new \code{lslxModel} object. ##
lslxModel$set("public",
              "initialize",
              function(model,
                       numeric_variable,
                       ordered_variable,
                       weight_variable,
                       auxiliary_variable,
                       group_variable,
                       reference_group,
                       level_group,
                       nlevel_ordered) {
                private$initialize_model(model = model)
                private$initialize_specification()
                private$initialize_tag(
                  numeric_variable = numeric_variable,
                  ordered_variable = ordered_variable,
                  weight_variable = weight_variable,
                  auxiliary_variable = auxiliary_variable,
                  group_variable = group_variable,
                  reference_group = reference_group,
                  level_group = level_group,
                  nlevel_ordered = nlevel_ordered
                )
                private$organize_specification()
                private$expand_specification()
              })

## \code{$initialize_model()} initializes a cleaned model. ##
lslxModel$set("private",
              "initialize_model",
              function(model) {
                model <-
                  gsub(pattern = "[[:blank:]]",
                       replacement = "",
                       x = model)
                model <-
                  gsub(pattern = "\\$|\\?|\\\\|\\^|%|&|#|\\[|\\]|\\{|\\}",
                       replacement = "",
                       x = model)
                model <-
                  gsub(pattern = ";",
                       replacement = "\n",
                       x = model)
                model <-
                  gsub(pattern = "\n{2,}",
                       replacement = "\n",
                       x = model)
                model <-
                  unlist(x = strsplit(x = model,
                                      split = "\n"),
                         use.names = FALSE)
                model <-
                  gsub(pattern = "=~",
                       replacement = ":=>",
                       x = model)
                model <-
                  gsub(pattern = "~~",
                       replacement = "<=>",
                       x = model)
                model <-
                  ifelse(
                    !grepl(pattern = "\\|=|=\\||\\|~|~\\|",
                           x = model),
                    gsub(
                      pattern = "\\|",
                      replacement = "|=",
                      x = model
                    ),
                    model
                  )
                model <-
                  ifelse(
                    !grepl(pattern = "<~|<~:|~>|:~>|<~>|\\|~|~\\|",
                           x = model),
                    gsub(
                      pattern = "~",
                      replacement = "<=",
                      x = model
                    ),
                    model
                  )
                self$model <- model
              })


## \code{$initialize_specification()} initializes a specification table. ##
lslxModel$set("private",
              "initialize_specification",
              function() {
                operator <- c("|=",
                              "=|",
                              "|~",
                              "~|",
                              "<=:",
                              "<~:",
                              ":=>",
                              ":~>",
                              "<=",
                              "<~",
                              "=>",
                              "~>",
                              "<=>",
                              "<~>")
                self$specification <-
                  do.call(what = rbind,
                          args = lapply(
                            X = self$model,
                            FUN = function(model_i) {
                              operator_i <- operator[sapply(
                                X = c(
                                  "\\|=",
                                  "=\\|",
                                  "\\|~",
                                  "~\\|",
                                  "<=:",
                                  "<~:",
                                  ":=>",
                                  ":~>",
                                  "<=[^:>]",
                                  "<~[^:>]",
                                  "[^:<]=>",
                                  "[^:<]~>",
                                  "<=>",
                                  "<~>"
                                ),
                                FUN = function(pattern) {
                                  grepl(pattern, model_i)
                                }
                              )]
                              if (length(operator_i) > 0) {
                                model_i_split <-
                                  strsplit(x = model_i, split = operator_i, fixed = TRUE)[[1]]
                                if (operator_i %in%
                                    c(":=>", ":~>", "=>", "~>", "=|", "~|")) {
                                  if (operator_i %in%
                                      c(":=>", ":~>", "=>", "~>")) {
                                    operator_i <-
                                      paste0(rev(gsub(
                                        pattern = ">",
                                        replacement = "<",
                                        x = substring(operator_i,
                                                      1:nchar(operator_i),
                                                      1:nchar(operator_i))
                                      )),
                                      collapse = "")
                                  } else {
                                    operator_i <-
                                      paste0(rev(substring(
                                        operator_i,
                                        1:nchar(operator_i),
                                        1:nchar(operator_i)
                                      )),
                                      collapse = "")
                                  }
                                  model_i <-
                                    c(left = model_i_split[2],
                                      operator = operator_i,
                                      right = model_i_split[1])
                                } else {
                                  model_i <-
                                    c(left = model_i_split[1],
                                      operator = operator_i,
                                      right = model_i_split[2])
                                }
                                left_i <-
                                  unlist(strsplit(x = model_i[["left"]],
                                                  split = "\\+"),
                                         use.names = FALSE)
                                right_i <-
                                  unlist(strsplit(x = model_i[["right"]],
                                                  split = "\\+"),
                                         use.names = FALSE)
                                model_i <-
                                  expand.grid(
                                    relation = NA_character_,
                                    left = left_i,
                                    right = right_i,
                                    operator =  operator_i,
                                    KEEP.OUT.ATTRS = FALSE,
                                    stringsAsFactors = FALSE
                                  )
                                left_i_split <-
                                  strsplit(model_i$left, split = "\\*")
                                right_i_split <-
                                  strsplit(model_i$right, split = "\\*")
                                model_i$left <-
                                  sapply(
                                    left_i_split,
                                    FUN = function(i)
                                      getElement(i, length(i))
                                  )
                                model_i$right <-
                                  sapply(
                                    right_i_split,
                                    FUN = function(i)
                                      getElement(i, length(i))
                                  )
                                if (any(model_i$right[model_i$operator %in%
                                                      c("<=>", "<~>")] == "1") |
                                    any(model_i$left[model_i$operator %in%
                                                     c("<=>", "<~>")] == "1")) {
                                  stop(
                                    "Intercept term '1' cannot present at any side of expression for covariance specification."
                                  )
                                }
                                if (any(model_i$left[!(model_i$operator %in%
                                                       c("<=>", "<~>"))] == "1")) {
                                  stop("Intercept term '1' cannot present at the arrow side of expression.")
                                }
                                model_i$left_prefix <-
                                  sapply(
                                    left_i_split,
                                    FUN = function(i) {
                                      ifelse(length(i) == 1L, NA_character_, i[1])
                                    }
                                  )
                                model_i$right_prefix <-
                                  sapply(
                                    right_i_split,
                                    FUN = function(i) {
                                      ifelse(length(i) == 1L, NA_character_, i[1])
                                    }
                                  )
                                if (any(!(is.na(model_i$left_prefix) |
                                          is.na(model_i$right_prefix)))) {
                                  stop("Prefix before '*' cannot simultaneously present at both side of expression.")
                                } else if (any(!is.na(model_i$left_prefix))) {
                                  model_i$prefix <-
                                    model_i$left_prefix
                                } else if (any(!is.na(model_i$right_prefix))) {
                                  model_i$prefix <-
                                    model_i$right_prefix
                                } else {
                                  model_i$prefix <- NA_character_
                                }
                                model_i$left_prefix <- NULL
                                model_i$right_prefix <- NULL
                                model_i$relation <-
                                  paste0(model_i$left,
                                         ifelse(
                                           model_i$operator %in%
                                             c("<=:", "<~:", "<=", "<~"),
                                           "<-",
                                           ifelse(model_i$operator %in%
                                                    c("<=>", "<~>"),
                                                  "<->",
                                                  "|")
                                         ),
                                         model_i$right)
                                model_i$operator <-
                                  ifelse(
                                    model_i$right == "1",
                                    ifelse(
                                      model_i$operator == "<=:",
                                      "<=",
                                      ifelse(model_i$operator == "<~:",
                                             "<~",
                                             model_i$operator)
                                    ),
                                    model_i$operator
                                  )
                              } else {
                                model_i <- data.frame()
                              }
                              return(model_i)
                            }
                          ))
                if (any(grepl(
                  pattern = "[[:digit:]]",
                  x = substr(
                    x = setdiff(x = unique(
                      c(self$specification$left,
                        self$specification$right)
                    ),
                    y = c("1")),
                    start = 1,
                    stop = 1
                  )
                ))) {
                  stop(
                    "Names of variable(s) or factor(s) cannot start with numbers.",
                    "\n  Please check the specified 'model'."
                  )
                }
              })


## \code{$initialize_tag()} initializes tags for variables. ##
lslxModel$set("private",
              "initialize_tag",
              function(numeric_variable,
                       ordered_variable,
                       weight_variable,
                       auxiliary_variable,
                       group_variable,
                       reference_group,
                       level_group,
                       nlevel_ordered) {
                self$numeric_variable <- numeric_variable
                self$ordered_variable <- ordered_variable
                self$weight_variable <- weight_variable
                self$auxiliary_variable <- auxiliary_variable
                self$group_variable <- group_variable
                self$reference_group <- reference_group
                self$level_group <- level_group
                self$name_factor <-
                  unique(self$specification[self$specification$operator %in%
                                              c("<=:", "<~:"),
                                            "right"])
                self$name_response <-
                  setdiff(x = unique(unlist(self$specification[!(self$specification$operator %in%
                                                                   c("|=", "|~")),
                                                               c("left", "right")])),
                          y = c(self$name_factor, "1"))
                
                if (!all(self$name_response %in% 
                         union(x = self$numeric_variable,
                               y = self$ordered_variable))) {
                  stop(
                    "Some response variable in 'model' cannot be found in 'data' or 'sample_cov'.",
                    "\n  Response variables specified by 'model' are ",
                    do.call(paste, as.list(self$name_response)),
                    ".",
                    "\n  Column names of 'data' or 'sample_cov' are ",
                    do.call(paste, as.list(union(x = self$numeric_variable,
                                                 y = self$ordered_variable))),
                    "."
                  )
                } 
                self$numeric_variable <- 
                  intersect(x = self$name_response,
                            y = self$numeric_variable)
                self$ordered_variable <- 
                  intersect(x = self$name_response,
                            y = self$ordered_variable)
                if (length(self$ordered_variable) > 0) {
                  self$nlevel_ordered <- 
                    nlevel_ordered[self$ordered_variable] 
                } else {
                  self$nlevel_ordered <- numeric(0)
                }
                self$auxiliary_variable <-
                  setdiff(x = self$auxiliary_variable,
                          y = self$name_response)
                self$name_eta <-
                  c(self$name_response, self$name_factor)
                self$name_endogenous <-
                  unique(self$specification[(self$specification$operator %in%
                                               c("<=:", "<~:", "<=", "<~")) &
                                              (self$specification$right != "1"),
                                            "left"])
                self$name_exogenous <-
                  unique(setdiff(x = self$name_eta,
                                 y = self$name_endogenous))
              })


## \code{$organize_specification()} organizes the specification table. ##
lslxModel$set("private",
              "organize_specification",
              function() {
                self$specification <-
                  do.call(what = rbind,
                          args = lapply(
                            split(self$specification, 1:nrow(self$specification)),
                            FUN = function(specification_i) {
                              if ((specification_i[["operator"]] %in% c("<=>", "<~>")) &
                                  (
                                    match(specification_i[["left"]], self$name_eta) <
                                    match(specification_i[["right"]], self$name_eta)
                                  )) {
                                specification_i <-
                                  data.frame(
                                    relation = paste0(specification_i[["right"]],
                                                      "<->",
                                                      specification_i[["left"]]),
                                    left = specification_i[["right"]],
                                    right = specification_i[["left"]],
                                    operator = specification_i[["operator"]],
                                    prefix = specification_i[["prefix"]]
                                  )
                              }
                              return(specification_i)
                            }
                          ))
                self$specification <-
                  self$specification[!duplicated(self$specification$relation,
                                                 fromLast = TRUE),]
                if (length(self$ordered_variable) > 0) {
                  relation_gamma <- 
                    setdiff(x = mapply(
                      FUN = function(ordered_variable_i, 
                                     nlevel_ordered_i) {
                        paste0(ordered_variable_i,
                               "|",
                               paste0("t", 1:(nlevel_ordered_i - 1)))
                    },
                    self$ordered_variable,
                    self$nlevel_ordered),
                    y = self$specification$relation)
                } else {
                  relation_gamma <- character(0)
                }
                if (length(relation_gamma) > 0) {
                  specification_gamma <-
                    data.frame(
                      relation = relation_gamma,
                      left = substr(
                        relation_gamma,
                        start = 1,
                        stop = regexpr("\\|", relation_gamma) - 1
                      ),
                      right = substr(
                        relation_gamma,
                        start = regexpr("\\|", relation_gamma) + 1,
                        stop = nchar(relation_gamma)
                      ),
                      operator = "|=",
                      prefix = NA_character_,
                      stringsAsFactors = FALSE
                    )
                } else {
                  specification_gamma = data.frame()
                }
                if (any(self$specification$right == "1")) {
                  if (length(intersect(x = self$numeric_variable,
                                       y = self$name_exogenous)) > 0) {
                    relation_alpha <-
                      setdiff(
                        x = paste(
                          intersect(
                            x = self$numeric_variable,
                            y = self$name_exogenous
                          ),
                          "1",
                          sep = "<-"
                        ),
                        y = self$specification$relation
                      )
                  } else {
                    relation_alpha <- character()
                  }
                } else {
                  if (length(self$numeric_variable) > 0) {
                    relation_alpha <-
                      setdiff(
                        x = paste(self$numeric_variable,
                                  "1",
                                  sep = "<-"),
                        y = self$specification$relation
                      )
                  } else {
                    relation_alpha <- character()
                  }
                }
                if (length(relation_alpha) > 0) {
                  specification_alpha <- data.frame(
                    relation = relation_alpha,
                    left = substr(
                      relation_alpha,
                      start = 1,
                      stop = regexpr("<-", relation_alpha) - 1
                    ),
                    right = "1",
                    operator = "<=",
                    prefix = NA_character_,
                    stringsAsFactors = FALSE
                  )
                } else {
                  specification_alpha = data.frame()
                }
                if (length(self$name_exogenous) > 1) {
                  relation_phi <-
                    setdiff(x = c(
                      paste0(self$name_eta,
                             "<->",
                             self$name_eta),
                      paste0(apply(
                        combn(rev(self$name_exogenous), 2),
                        2,
                        FUN = function(i) {
                          paste(i[1], i[2], sep = "<->")
                        }
                      ))
                    ),
                    y = self$specification$relation)
                } else {
                  relation_phi <-
                    setdiff(
                      x = paste0(self$name_eta,
                                 "<->",
                                 self$name_eta),
                      y = self$specification$relation
                    )
                }
                if (length(relation_phi) > 0) {
                  specification_phi <-
                    data.frame(
                      relation = relation_phi,
                      left = substr(
                        relation_phi,
                        start = 1,
                        stop = regexpr("<->", relation_phi) - 1
                      ),
                      right = substr(
                        relation_phi,
                        start = regexpr("<->", relation_phi) + 3,
                        stop = nchar(relation_phi)
                      ),
                      operator = "<=>",
                      prefix = NA_character_,
                      stringsAsFactors = FALSE
                    )
                } else {
                  specification_phi = data.frame()
                }
                self$specification <-
                  rbind(
                    self$specification,
                    specification_gamma,
                    specification_alpha,
                    specification_phi,
                    make.row.names = FALSE,
                    stringsAsFactors = FALSE
                  )
              })



## \code{$expand_specification()} expands the specification table. ##
lslxModel$set("private",
              "expand_specification",
              function() {
                prefix_split <-
                  lapply(X = self$specification$prefix,
                         FUN = function(prefix_i) {
                           if (!is.na(prefix_i)) {
                             if ((substr(prefix_i, 
                                         start = 1, 
                                         stop = 2) == "c(") & 
                                 (substr(prefix_i, 
                                         start = nchar(prefix_i), 
                                         stop = nchar(prefix_i)) == ")")) {
                               prefix_i <- 
                                 substr(prefix_i, start = 3, stop = (nchar(prefix_i) - 1))
                               prefix_i <- strsplit(x = prefix_i, split = ",")[[1]]
                             }
                             prefix_i <- 
                               ifelse(gsub(pattern = "\\(.*$", replacement = "", x = prefix_i) %in% 
                                        c("free", "fix", "pen", "start", "lab"),
                                      prefix_i,
                                      ifelse(suppressWarnings(!is.na(as.numeric(prefix_i))),
                                             paste0("fix", "(",  prefix_i, ")"),
                                             paste0("lab", "(",  prefix_i, ")")))
                           } 
                           return(prefix_i)
                         })
                
                
                if (anyNA(self$reference_group)) {
                  if (any(sapply(X = prefix_split,
                                 FUN = function(prefix_split_i) {
                                   return(length(prefix_split_i) > 1)
                                 }))) {
                    stop("Vectorized prefix cannot be applied to the case of 'reference_group = NA'.")
                  }
                  self$specification <- 
                    do.call(what = rbind,
                            args = lapply(
                              X = c("<NA>", self$level_group),
                              FUN = function(level_group_i) {
                                specification_i <- self$specification
                                specification_i$group <- level_group_i
                                specification_i$reference <-
                                  ifelse(level_group_i == "<NA>", TRUE, FALSE)
                                specification_i$operator <-
                                  ifelse(specification_i$group == "<NA>", 
                                         specification_i$operator, 
                                         gsub(pattern = "=",
                                              replacement = "~",
                                              x = specification_i$operator))
                                specification_i$prefix <-
                                  sapply(
                                    X = prefix_split,
                                    FUN = function(prefix_split_j) {
                                      if (level_group_i != "<NA>") {
                                        if (gsub(
                                          pattern = "\\(.*$", 
                                          replacement = "", 
                                          x = prefix_split_j) %in% "fix") {
                                          prefix_split_j <- "fix(0)"
                                        } else {
                                          prefix_split_j <- NA
                                        }
                                      } 
                                      return(prefix_split_j)
                                    }
                                  )
                                rownames(specification_i) <-
                                  paste0(specification_i$relation,
                                         "/",
                                         specification_i$group)
                                return(specification_i)
                              }
                            ))
                } else {
                  if (any(sapply(X = prefix_split,
                                 FUN = function(prefix_split_i) {
                                   return((length(prefix_split_i) > 1) & 
                                          (length(prefix_split_i) != length(self$level_group)))
                                 }))) {
                    stop("The length of prefix vector should be 1 or equal to to the number of groups.")
                  }
                  self$specification <- 
                    do.call(what = rbind,
                            args = lapply(
                              X = self$level_group,
                              FUN = function(level_group_i) {
                                specification_i <- self$specification
                                specification_i$group <- level_group_i
                                specification_i$reference <-
                                  ifelse(
                                    is.null(self$reference_group),
                                    FALSE,
                                    ifelse(level_group_i == self$reference_group,
                                           TRUE,
                                           FALSE)
                                  )
                                specification_i$prefix <-
                                  sapply(
                                    X = prefix_split,
                                    FUN = function(prefix_split_j) {
                                      if (length(prefix_split_j) == 1) {
                                        if (!is.null(self$reference_group)) {
                                          if (level_group_i != self$reference_group) {
                                            prefix_split_j_verb <- 
                                              gsub(pattern = "\\(.*$", 
                                                   replacement = "", 
                                                   x = prefix_split_j)
                                            if (prefix_split_j_verb %in% 
                                                c("free", "fix", "pen", "start")) {
                                              prefix_split_j <- 
                                                paste0(prefix_split_j_verb, "(", 0, ")")
                                            } else {
                                              prefix_split_j <- NA
                                            }
                                          } 
                                        }
                                      } else {
                                        prefix_split_j <-
                                          prefix_split_j[level_group_i == self$level_group]
                                      }
                                      return(prefix_split_j)
                                    }
                                  )
                                rownames(specification_i) <-
                                  paste0(specification_i$relation,
                                         "/",
                                         specification_i$group)
                                return(specification_i)
                              }
                            ))
                }
                self$specification$matrix <-
                  factor(
                    x = ifelse(
                      self$specification$operator %in% c("|=", "|~"),
                      "gamma",
                      ifelse(
                        self$specification$operator %in% c("<=>", "<~>"),
                        "phi",
                        ifelse(self$specification$right == "1",
                               "alpha",
                               "beta")
                      )
                    ),
                    levels = c("gamma", "alpha", "beta", "phi")
                  )
                self$specification$block <-
                  with(self$specification, {
                    block_left <-
                      ifelse(left %in% self$name_response,
                             "y",
                             "f")
                    block_right <-
                      ifelse(matrix %in% "gamma",
                             "t",
                             ifelse(
                               right %in% self$name_response,
                               "y",
                               ifelse(
                                 right %in% self$name_factor,
                                 "f",
                                 "1"
                               )
                             ))
                    block_middle <-
                      ifelse(matrix %in% "gamma",
                             "|",
                             ifelse(matrix %in% c("alpha", "beta"),
                                    "<-",
                                    "<->"))
                    paste0(block_left, block_middle, block_right)
                  })
                
                self$specification$type <-
                  with(self$specification, {
                    prefix_verb <- 
                      gsub(pattern = "\\(.*$", replacement = "", x = prefix)
                    type <- 
                      ifelse(prefix_verb %in% "free",
                             "free",
                             ifelse(prefix_verb %in% "fix",
                                    "fixed",
                                    ifelse(prefix_verb %in% "pen",
                                           "pen",
                                           ifelse(
                                             operator %in% c("|=", "<=:", "<=", "<=>"),
                                             "free",
                                             "pen"))))
                    return(type)
                  })
                self$specification$start <-
                  with(self$specification, {
                    prefix_verb <- 
                      gsub(pattern = "\\(.*$", replacement = "", x = prefix)
                    prefix_value <- 
                      substr(x = prefix, 
                             start = (regexpr("\\(", prefix) + 1), 
                             stop = (nchar(prefix) - 1))
                    start <- 
                      ifelse(prefix_verb %in% c("free", "fix", "pen", "start"),
                             suppressWarnings(as.numeric(prefix_value)),
                             NA_real_)
                    return(start)
                  })
                
                self$specification$label <-
                  with(self$specification, {
                    prefix_verb <- 
                      gsub(pattern = "\\(.*$", replacement = "", x = prefix)
                    prefix_value <- 
                      substr(x = prefix, 
                             start = (regexpr("\\(", prefix) + 1), 
                             stop = (nchar(prefix) - 1))
                    label <- 
                      ifelse(prefix_verb %in% "lab",
                             prefix_value,
                             NA_character_)
                    return(label)
                  })
                self$specification$operator <- NULL
                self$specification$prefix <- NULL
                self$specification <-
                  self$specification[order(
                    self$specification$reference,
                    self$specification$group,
                    self$specification$matrix,
                    self$specification$block,
                    match(self$specification$right, self$name_eta),
                    match(self$specification$left, self$name_eta),
                    method = "radix"
                  ),]
                self$specification <-
                  self$specification[order(
                    self$specification$reference,
                    decreasing = TRUE,
                    method = "radix"
                  ), ]
                
              })
