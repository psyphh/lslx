## define R6 class \code{lslxModel} to store model information. ##
lslxModel <-
  R6::R6Class(
    classname = "lslxModel",
    public = list(
      ordered_variable = "character",
      weight_variable = "character",
      auxiliary_variable = "character",
      group_variable = "character",
      reference_group = "character",
      name_group = "character",
      name_response = "character",
      name_factor = "character",
      name_eta = "character",
      name_exogenous = "character",
      name_endogenous = "character",
      specification = "data.frame"
    )
  )

## \code{$new()} initializes a new \code{lslxModel} object. ##
lslxModel$set("public",
              "initialize",
              function(model,
                       ordered_variable,
                       weight_variable,
                       auxiliary_variable,
                       group_variable,
                       reference_group,
                       name_group) {
                private$initialize_specification(model = model)
                private$initialize_variable(ordered_variable = ordered_variable,
                                            weight_variable = weight_variable,
                                            auxiliary_variable = auxiliary_variable,
                                            group_variable = group_variable,
                                            reference_group = reference_group,
                                            name_group = name_group)
                private$expand_specification()
                private$expand_specification_alpha()
                private$expand_specification_phi()
                self$specification <-
                  self$specification[order(
                    self$specification$reference,
                    self$specification$group,
                    self$specification$matrix,
                    self$specification$block,
                    match(self$specification$right, self$name_eta),
                    match(self$specification$left, self$name_eta),
                    decreasing = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
                    method = "radix"
                  ),]
              })

## \code{$initialize_specification()} initializes a specification table. ##
lslxModel$set("private",
              "initialize_specification",
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
                    gsub(pattern = "\\|",
                         replacement = "|=",
                         x = model),
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
                model <-
                  sapply(
                    X = c("\\|=", "=\\|", "\\|~", "~\\|",
                          "<=:", "<~:", ":=>", ":~>",
                          "<=", "<~", "=>", "~>", 
                          "<=>", "<~>"),
                    FUN = function(operator_i) {
                      idx <-
                        grep(
                          pattern = ifelse(
                            operator_i %in% c("=>", "~>"),
                            paste0("[^:<]", operator_i),
                            ifelse(operator_i %in% c("<=", "<~"),
                                   paste0(operator_i, "[^:>]"),
                                   operator_i)
                          ),
                          x = model,
                          value = FALSE
                        )
                      model_i <-
                        do.call(what = rbind,
                                args = strsplit(x = model[idx],
                                                split = operator_i))
                      if (length(model_i) > 0) {
                        if (operator_i %in% c(":=>", ":~>", "=>", "~>", "=\\|", "~\\|")) {
                          model_i <-
                            model_i[, c(2, 1), drop = FALSE]
                          if (operator_i %in% c(":=>", ":~>", "=>", "~>")) {
                            operator_i <-
                              paste0(rev(gsub(
                                pattern = ">",
                                replacement = "<",
                                x = substring(operator_i,
                                              1:nchar(operator_i),
                                              1:nchar(operator_i))
                              )),
                              collapse = "")
                          } else if (operator_i %in% c("=\\|", "~\\|")) {
                            operator_i <-
                            paste0(rev(substring(operator_i,
                                                 1:nchar(operator_i),
                                                 1:nchar(operator_i))),
                            collapse = "")
                          } else {
                            
                          }
                        }
                        model_i <-
                          cbind(idx = idx,
                                "colnames<-"(model_i,
                                             value = c("left", "right")),
                                operator = gsub("\\\\", "", operator_i))
                      }
                      return(model_i)
                    },
                    simplify = TRUE,
                    USE.NAMES = TRUE
                  )
                model <-
                  do.call(what = rbind, args = model)
                model <-
                  model[order(as.numeric(model[, "idx"])),
                              ,
                              drop = FALSE]
                model <- model[, -1, drop = FALSE]
                self$specification <-
                  do.call(what = rbind, 
                          args = apply(
                            X = model,
                            MARGIN = 1,
                            FUN = function(model_i) {
                              model_i_left <-
                                unlist(strsplit(x = model_i[["left"]],
                                                split = "\\+"),
                                       use.names = FALSE)
                              model_i_right <-
                                unlist(strsplit(x = model_i[["right"]],
                                                split = "\\+"),
                                       use.names = FALSE)
                              model_i <-
                                expand.grid(
                                  relation = NA_character_,
                                  left = model_i_left,
                                  right = model_i_right,
                                  operator = model_i["operator"],
                                  KEEP.OUT.ATTRS = FALSE,
                                  stringsAsFactors = FALSE
                                )
                              model_i_left_split <-
                                strsplit(model_i$left, split = "\\*")
                              model_i_right_split <-
                                strsplit(model_i$right, split = "\\*")
                              model_i$left <-
                                sapply(
                                  model_i_left_split,
                                  FUN = function(i)
                                    getElement(i, length(i))
                                )
                              model_i$right <-
                                sapply(
                                  model_i_right_split,
                                  FUN = function(i)
                                    getElement(i, length(i))
                                )
                              if (any(model_i$right[
                                model_i$operator %in% c("<=>", "<~>")] == "1") | 
                                any(model_i$left[
                                  model_i$operator %in% c("<=>", "<~>")] == "1")) {
                                stop(
                                  "Intercept term '1' cannot present at the expression for covariance.",
                                  "\n  Please check the specified 'model'."
                                )
                              }
                              if (any(model_i$left[
                                !(model_i$operator %in% c("<=>", "<~>"))] == "1")) {
                                stop(
                                  "Intercept term '1' cannot present at the arrow side of expression.",
                                  "\n  Please check the specified 'model'."
                                )
                              }
                              model_i$left_prefix <-
                                sapply(
                                  model_i_left_split,
                                  FUN = function(i) {
                                    ifelse(length(i) == 1L, NA_character_, i[1])
                                  }
                                )
                              model_i$right_prefix <-
                                sapply(
                                  model_i_right_split,
                                  FUN = function(i) {
                                    ifelse(length(i) == 1L, NA_character_, i[1])
                                  }
                                )
                              if (any(!(is.na(model_i$left_prefix) |
                                        is.na(model_i$right_prefix)))) {
                                stop(
                                  "Prefix before '*' cannot simultaneously present at both side of expression.",
                                  "\n  Please check the specified 'model'."
                                )
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
                                paste0(
                                  model_i$left,
                                  ifelse(
                                    model_i$operator %in% 
                                      c("<=:", "<~:", "<=", "<~"),
                                    "<-",
                                    ifelse(model_i$operator %in% 
                                             c("<=>", "<~>"),
                                           "<->",
                                           "|")),
                                  model_i$right
                                )
                              model_i$operator <-
                                ifelse(
                                  model_i$right == "1",
                                  ifelse(
                                    model_i$operator == "<=:",
                                    "<=",
                                    ifelse(
                                      model_i$operator == "<~:",
                                      "<~",
                                      model_i$operator
                                    )
                                  ),
                                  model_i$operator
                                )
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


## \code{$initialize_variable()} initializes variable names. ##
lslxModel$set("private",
              "initialize_variable",
              function(ordered_variable,
                       weight_variable,
                       auxiliary_variable,
                       group_variable,
                       reference_group,
                       name_group) {
                self$ordered_variable <- ordered_variable
                self$weight_variable <- weight_variable
                self$auxiliary_variable <- auxiliary_variable
                self$group_variable <- group_variable
                self$reference_group <- reference_group
                self$name_group <- name_group
                self$name_factor <-
                  unique(self$specification[
                    self$specification$operator %in% c("<=:", "<~:"), "right"])
                self$name_response <-
                  setdiff(
                    x = unique(unlist(self$specification[
                      !(self$specification$operator %in% c("|=", "|~")), 
                      c("left", "right")])),
                    y = c(self$name_factor, "1"))
                self$auxiliary_variable <-
                  setdiff(x = self$auxiliary_variable,
                          y = self$name_response)
                self$name_eta <-
                  c(self$name_response, self$name_factor)
                self$name_endogenous <-
                  unique(self$specification[
                    (self$specification$operator %in% c("<=:", "<~:", "<=", "<~")) &
                      (self$specification$right != "1"), "left"])
                self$name_exogenous <-
                  unique(setdiff(x = self$name_eta,
                                 y = self$name_endogenous))
              })


## \code{$expand_specification()} expands a specification table. ##
lslxModel$set("private",
              "expand_specification",
              function() {
                self$specification <-
                  data.frame(t(apply(
                    self$specification,
                    1,
                    FUN = function(specification_i) {
                      if ((specification_i[["operator"]] %in% c("<=>", "<~>")) &
                          (match(specification_i[["left"]], self$name_eta) <
                           match(specification_i[["right"]], self$name_eta))) {
                        specification_i <-
                          c(
                            relation = paste0(specification_i[["right"]],
                                              "<->",
                                              specification_i[["left"]]),
                            left = specification_i[["right"]],
                            operator = specification_i[["operator"]],
                            right = specification_i[["left"]],
                            prefix = specification_i[["prefix"]]
                          )
                      } else {
                        specification_i <-
                          c(
                            relation = specification_i[["relation"]],
                            left = specification_i[["left"]],
                            operator = specification_i[["operator"]],
                            right = specification_i[["right"]],
                            prefix = specification_i[["prefix"]]
                          )
                      }
                      return(specification_i)
                    }
                  )), 
                  stringsAsFactors = FALSE)
                
                self$specification <-
                  self$specification[!duplicated(self$specification$relation,
                                                 fromLast = TRUE),]
                self$specification$prefix <- 
                  gsub(pattern = " ", replacement = "",
                       x = self$specification$prefix)
                prefix_split <-
                  strsplit(self$specification$prefix, ",")
                self$specification <-
                  do.call(what = rbind,
                          args = lapply(
                            X = self$name_group,
                            FUN = function(name_group_i) {
                              specification_i <- self$specification
                              specification_i$group <- name_group_i
                              specification_i$reference <-
                                ifelse(
                                  is.null(self$reference_group),
                                  FALSE,
                                  ifelse(name_group_i == self$reference_group,
                                         TRUE,
                                         FALSE)
                                )
                              specification_i$prefix <-
                                sapply(
                                  X = prefix_split,
                                  FUN = function(prefix_split_j) {
                                    if (length(prefix_split_j) == 1) {
                                      if (!is.null(self$reference_group)) {
                                        if (name_group_i == self$reference_group) {
                                          prefix_split_j <- prefix_split_j
                                        } else {
                                          prefix_split_j <-
                                            gsub(pattern = "[[:digit:].-]",
                                                 replacement = 0,
                                                 x = prefix_split_j)
                                        }
                                      }
                                      prefix_split_j <-
                                        prefix_split_j
                                    } else {
                                      if (length(prefix_split_j) ==
                                          length(self$name_group)) {
                                        prefix_split_j <-
                                          prefix_split_j[name_group_i == self$name_group]
                                      } else {
                                        stop("The length of prefix vector should be 1 or equal to to the number of group(s)")
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
                self$specification$matrix <-
                  ifelse(self$specification$operator %in% c("|=", "|~"),
                         "gamma",
                         ifelse(self$specification$operator %in% c("<=>", "<~>"),
                                "phi",
                                ifelse(self$specification$right == "1",
                                       "alpha",
                                       "beta")))
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
                  ifelse(
                    grepl("[[:digit:]]", self$specification$prefix) &
                      !grepl("[[:alpha:]]", self$specification$prefix),
                    "fixed",
                    ifelse(
                      grepl("fix", self$specification$prefix),
                      "fixed",
                      ifelse(
                        grepl("pen", self$specification$prefix),
                        "pen",
                        ifelse(
                          grepl("free", self$specification$prefix),
                          "free",
                          ifelse(
                            self$specification$operator %in%
                              c("<=:", "<=", "<=>"),
                            "free",
                            "pen"
                          )
                        )
                      )
                    )
                  )
                self$specification$start <-
                  as.numeric(gsub("[^[:digit:].-]",
                                  "",
                                  self$specification$prefix))
                if (any(grepl("fix", self$specification$prefix) &
                        is.na(self$specification$start))) {
                  stop("Some fixed coefficient contains 'NA' or unspecified starting value.")
                }
                self$specification$operator <- NULL
                self$specification$prefix <- NULL
              })

## \code{$expand_specification_alpha()} expand the specification table for intercepts. ##
lslxModel$set("private",
              "expand_specification_alpha",
              function() {
                specification_alpha <-
                  do.call(what = rbind.data.frame,
                          args = lapply(
                            X = self$name_group,
                            FUN = function(name_group_i) {
                              if (any(self$specification$matrix[
                                self$specification$group == name_group_i] == "alpha")) {
                                if (length(intersect(
                                  x = self$name_response,
                                  y = self$name_exogenous
                                )) > 0) {
                                  relation_alpha_i <-
                                    setdiff(
                                      x = paste(
                                        intersect(
                                          x = self$name_response,
                                          y = self$name_exogenous
                                        ),
                                        "1",
                                        sep = "<-"
                                      ),
                                      y = self$specification$relation[self$specification$group == name_group_i]
                                    )
                                } else {
                                  relation_alpha_i <- character()
                                }
                              } else {
                                relation_alpha_i <-
                                  setdiff(
                                    x = paste(self$name_response,
                                              "1",
                                              sep = "<-"),
                                    y = self$specification$relation[self$specification$group == name_group_i]
                                  )
                              }
                              if (length(relation_alpha_i) > 0) {
                                specification_alpha_i <- data.frame(
                                  relation = relation_alpha_i,
                                  left = substr(
                                    relation_alpha_i,
                                    start = 1,
                                    stop = regexpr("<-", relation_alpha_i) - 1
                                  ),
                                  right = "1",
                                  group = name_group_i,
                                  reference = ifelse(
                                    is.null(self$reference_group),
                                    FALSE,
                                    ifelse(name_group_i == self$reference_group,
                                           TRUE,
                                           FALSE)
                                  ),
                                  matrix = "alpha",
                                  block = "y<-1",
                                  type = "free",
                                  start = NA_real_,
                                  stringsAsFactors = FALSE
                                )
                                rownames(specification_alpha_i) <-
                                  paste0(specification_alpha_i$relation,
                                         "/",
                                         specification_alpha_i$group)
                              } else {
                                specification_alpha_i = data.frame()
                              }
                              return(specification_alpha_i)
                            }
                          ))
                self$specification <-
                  rbind(self$specification,
                        specification_alpha,
                        stringsAsFactors = FALSE)
              })

## \code{$expand_specification_phi()} expand the specification table for covariances. ##
lslxModel$set("private",
              "expand_specification_phi",
              function() {
                specification_phi <-
                  do.call(what = rbind.data.frame,
                          args = lapply(
                            X = self$name_group,
                            FUN = function(name_group_i) {
                              if (length(self$name_exogenous) > 1) {
                                relation_phi_i <-
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
                                  y = self$specification$relation[self$specification$group == name_group_i])
                              } else {
                                relation_phi_i <-
                                  setdiff(
                                    x = paste0(self$name_eta,
                                               "<->",
                                               self$name_eta),
                                    y = self$specification$relation[self$specification$group == name_group_i]
                                  )
                              }
                              if (length(relation_phi_i) > 1) {
                                left_phi_i <-
                                  substr(relation_phi_i,
                                         start = 1,
                                         stop = regexpr("<->", relation_phi_i) - 1)
                                right_phi_i <-
                                  substr(
                                    relation_phi_i,
                                    start = regexpr("<->", relation_phi_i) + 3,
                                    stop = nchar(relation_phi_i)
                                  )
                                specification_phi_i <-
                                  data.frame(
                                    relation = relation_phi_i,
                                    left = left_phi_i,
                                    right = right_phi_i,
                                    group = name_group_i,
                                    reference = ifelse(
                                      is.null(self$reference_group),
                                      FALSE,
                                      ifelse(name_group_i == self$reference_group,
                                             TRUE,
                                             FALSE)
                                    ),
                                    matrix = "phi",
                                    block =
                                      ifelse(
                                        left_phi_i %in% self$name_factor,
                                        ifelse(right_phi_i %in% self$name_factor,
                                               "f<->f",
                                               "f<->y"),
                                        ifelse(right_phi_i %in% self$name_factor,
                                               "y<->f",
                                               "y<->y")
                                      ),
                                    type = "free",
                                    start = NA_real_,
                                    stringsAsFactors = FALSE
                                  )
                                rownames(specification_phi_i) <-
                                  paste(specification_phi_i$relation,
                                        specification_phi_i$group,
                                        sep = "/")
                                return(specification_phi_i)
                              }
                            }
                          ))
                self$specification <-
                  rbind(self$specification,
                        specification_phi,
                        stringsAsFactors = FALSE)
              })
