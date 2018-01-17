## define R6 class \code{lslxModel} to store model information. ##
lslxModel <-
  R6::R6Class(
    classname = "lslxModel",
    public = list(
      name_group = "character",
      reference_group = "character",
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
                       name_group,
                       reference_group) {
                self$name_group <- name_group
                self$reference_group <- reference_group
                model_parsed <-
                  private$parse_model(model = model)
                self$name_factor <-
                  unique(model_parsed[model_parsed$operator %in%
                                        c("<=:", "<~:"), "right"])
                self$name_response <-
                  setdiff(x = unique(unlist(model_parsed[, c("left", "right")])),
                          y = c(self$name_factor, "1"))
                self$name_eta <-
                  c(self$name_response, self$name_factor)
                self$name_endogenous <-
                  unique(model_parsed[(model_parsed$operator %in%
                                         c("<=:", "<~:", "<=", "<~")) &
                                        (model_parsed$right != "1"),
                                      "left"])
                self$name_exogenous <-
                  unique(setdiff(x = self$name_eta,
                                 y = self$name_endogenous))
                private$initialize_specification(model_parsed = model_parsed)
                private$expand_specification_alpha()
                private$expand_specification_psi()
                self$specification <-
                  self$specification[order(
                    self$specification$reference,
                    self$specification$group,
                    self$specification$matrice,
                    self$specification$block,
                    match(self$specification$right, self$name_eta),
                    match(self$specification$left, self$name_eta),
                    decreasing = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
                    method = "radix"
                  ), ]
              })

## \code{$parse_model()} parses specified model. ##
lslxModel$set("private",
              "parse_model",
              function(model) {
                model_cleaned <-
                  gsub(pattern = "[[:blank:]]",
                       replacement = "",
                       x = model)
                model_cleaned <-
                  gsub(pattern = "\\$|\\?|\\||\\\\|\\^|%|&|#|\\[|\\]|\\{|\\}",
                       replacement = "",
                       x = model_cleaned)
                model_cleaned <-
                  gsub(pattern = ";",
                       replacement = "\n",
                       x = model_cleaned)
                model_cleaned <-
                  gsub(pattern = "\n{2,}",
                       replacement = "\n",
                       x = model_cleaned)
                model_split <-
                  unlist(x = strsplit(x = model_cleaned,
                                      split = "\n"),
                         use.names = FALSE)
                model_split <-
                  sapply(
                    X = c("<=:", "<~:", ":=>", ":~>",
                          "<=", "<~", "=>", "~>", "<=>", "<~>"),
                    FUN = function(operator_i) {
                      idx <-
                        grep(
                          pattern = ifelse(
                            operator_i %in% c("=>", "~>"),
                            paste0("[^:<]", operator_i),
                            paste0(operator_i, "[^:>]")
                          ),
                          x = model_split,
                          value = FALSE
                        )
                      model_split_i <-
                        do.call(what = rbind,
                                args = strsplit(x = model_split[idx],
                                                split = operator_i))
                      if (length(model_split_i) > 0) {
                        if (operator_i %in% c(":=>", ":~>", "=>", "~>")) {
                          model_split_i <-
                            model_split_i[, c(2, 1), drop = FALSE]
                          operator_i <-
                            paste0(rev(gsub(
                              pattern = ">",
                              replacement = "<",
                              x = substring(operator_i,
                                            1:nchar(operator_i),
                                            1:nchar(operator_i))
                            )),
                            collapse = "")
                        }
                        model_split_i <-
                          cbind(idx = idx,
                                "colnames<-"(model_split_i,
                                             value = c("left", "right")),
                                operator = operator_i)
                      }
                      return(model_split_i)
                    },
                    simplify = TRUE,
                    USE.NAMES = TRUE
                  )
                model_split <-
                  do.call(what = rbind, args = model_split)
                model_split <-
                  model_split[order(as.numeric(model_split[, "idx"])),
                              ,
                              drop = FALSE]
                model_split <-
                  model_split[, -1, drop = FALSE]
                model_parsed <-
                  apply(
                    X = model_split,
                    MARGIN = 1,
                    FUN = function(model_split_i) {
                      left_split <-
                        unlist(strsplit(x = model_split_i[["left"]],
                                        split = "\\+"),
                               use.names = FALSE)
                      right_split <-
                        unlist(strsplit(x = model_split_i[["right"]],
                                        split = "\\+"),
                               use.names = FALSE)
                      model_parsed_i <-
                        expand.grid(
                          relation = NA_character_,
                          left = left_split,
                          right = right_split,
                          operator = model_split_i["operator"],
                          KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE
                        )
                      left_split <-
                        strsplit(model_parsed_i$left, split = "\\*")
                      right_split <-
                        strsplit(model_parsed_i$right, split = "\\*")
                      model_parsed_i$left <-
                        sapply(
                          left_split,
                          FUN = function(i)
                            getElement(i, length(i))
                        )
                      model_parsed_i$right <-
                        sapply(
                          right_split,
                          FUN = function(i)
                            getElement(i, length(i))
                        )
                      
                      #check invalid intercept specification
                      if ("|"(any(
                        "["(
                          model_parsed_i$right,
                          model_parsed_i$operator %in%
                          c("<=>", "<~>")
                        ) == "1"
                      ),
                      any(
                        "["(
                          model_parsed_i$left,
                          model_parsed_i$operator %in%
                          c("<=>", "<~>")
                        ) == "1"
                      ))) {
                        stop(
                          "Intercept term '1' cannot present at the expression for covariance.",
                          "\n  Please check the specified 'model'."
                        )
                      }
                      if (any("["(model_parsed_i$left, !(
                        model_parsed_i$operator %in% c("<=>", "<~>")
                      )) == "1")) {
                        stop(
                          "Intercept term '1' cannot present at the arrow side of expression.",
                          "\n  Please check the specified 'model'."
                        )
                      }
                      model_parsed_i$left_prefix <-
                        sapply(
                          left_split,
                          FUN = function(i) {
                            ifelse(length(i) == 1L, NA_character_, i[1])
                          }
                        )
                      model_parsed_i$right_prefix <-
                        sapply(
                          right_split,
                          FUN = function(i) {
                            ifelse(length(i) == 1L, NA_character_, i[1])
                          }
                        )
                      if (any(!(
                        is.na(model_parsed_i$left_prefix) |
                        is.na(model_parsed_i$right_prefix)
                      ))) {
                        stop(
                          "Prefix before '*' cannot simultaneously present at both side of expression.",
                          "\n  Please check the specified 'model'."
                        )
                      } else if (any(!is.na(model_parsed_i$left_prefix))) {
                        model_parsed_i$prefix <-
                          model_parsed_i$left_prefix
                      } else if (any(!is.na(model_parsed_i$right_prefix))) {
                        model_parsed_i$prefix <-
                          model_parsed_i$right_prefix
                      } else {
                        model_parsed_i$prefix <- NA_character_
                      }
                      model_parsed_i$left_prefix <- NULL
                      model_parsed_i$right_prefix <- NULL
                      model_parsed_i$relation <-
                        paste0(
                          model_parsed_i$left,
                          ifelse(
                            model_parsed_i$operator %in%
                              c("<=:", "<~:", "<=", "<~"),
                            "<-",
                            "<->"
                          ),
                          model_parsed_i$right
                        )
                      model_parsed_i$operator <-
                        ifelse(
                          model_parsed_i$right == "1",
                          ifelse(
                            model_parsed_i$operator == "<=:",
                            "<=",
                            ifelse(
                              model_parsed_i$operator == "<~:",
                              "<~",
                              model_parsed_i$operator
                            )
                          ),
                          model_parsed_i$operator
                        )
                      return(model_parsed_i)
                    }
                  )
                model_parsed <-
                  do.call(what = rbind, args = model_parsed)
                if (any(grepl(
                  pattern = "[[:digit:]]",
                  x = substr(
                    x = setdiff(x = unique(
                      c(model_parsed$left,
                        model_parsed$right)
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
                return(model_parsed)
              })

## \code{$initialize_specification()} initializes a specification table. ##
lslxModel$set("private",
              "initialize_specification",
              function(model_parsed) {
                model_parsed <-
                  apply(
                    model_parsed,
                    1,
                    FUN = function(model_parsed_i) {
                      if ((model_parsed_i[["operator"]] %in% c("<=>", "<~>")) &
                          (match(model_parsed_i[["left"]], self$name_eta) <
                           match(model_parsed_i[["right"]], self$name_eta))) {
                        model_parsed_i <-
                          c(
                            relation = paste0(model_parsed_i[["right"]],
                                              "<->",
                                              model_parsed_i[["left"]]),
                            left = model_parsed_i[["right"]],
                            operator = model_parsed_i[["operator"]],
                            right = model_parsed_i[["left"]],
                            prefix = model_parsed_i[["prefix"]]
                          )
                      } else {
                        model_parsed_i <-
                          c(
                            relation = model_parsed_i[["relation"]],
                            left = model_parsed_i[["left"]],
                            operator = model_parsed_i[["operator"]],
                            right = model_parsed_i[["right"]],
                            prefix = model_parsed_i[["prefix"]]
                          )
                      }
                      return(model_parsed_i)
                    }
                  )
                self$specification <-
                  data.frame(t(model_parsed), stringsAsFactors = FALSE)
                
                self$specification <-
                  self$specification[!duplicated(self$specification$relation,
                                                 fromLast = TRUE), ]
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
                                  is.na(self$reference_group),
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
                                      if (!is.na(self$reference_group)) {
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
                                       "|",
                                       specification_i$group)
                              return(specification_i)
                            }
                          ))
                self$specification$matrice <-
                  ifelse(self$specification$operator %in% c("<=>", "<~>"),
                         "psi",
                         ifelse(
                           !(self$specification$operator %in%
                               c("<=:", "<~:", "<=", "<~")),
                           "tau",
                           ifelse(self$specification$right == "1",
                                  "alpha",
                                  "beta")
                         ))
                self$specification$block <-
                  mapply(
                    self$specification$left,
                    self$specification$right,
                    self$specification$matrice,
                    FUN = function(left, right, matrice) {
                      block_left <-
                        ifelse(left %in% self$name_response,
                               "y",
                               "f")
                      block_right <-
                        ifelse(
                          right %in% self$name_response,
                          "y",
                          ifelse(
                            right %in% self$name_factor,
                            "f",
                            ifelse(right == "1",
                                   "1",
                                   "()")
                          )
                        )
                      block_middle <-
                        ifelse(matrice %in% c("alpha", "beta"),
                               "<-",
                               ifelse(matrice == "psi",
                                      "<->",
                                      ""))
                      block <-
                        paste0(block_left, block_middle, block_right)
                      return(block)
                    }
                  )
                self$specification$type <-
                  ifelse(grepl("fix", self$specification$prefix),
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
                         ))
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
                              if (any(self$specification$matrice[self$specification$group == name_group_i] == "alpha")) {
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
                                    is.na(self$reference_group),
                                    FALSE,
                                    ifelse(name_group_i == self$reference_group,
                                           TRUE,
                                           FALSE)
                                  ),
                                  matrice = "alpha",
                                  block = "y<-1",
                                  type = "free",
                                  start = NA_real_,
                                  stringsAsFactors = FALSE
                                )
                                rownames(specification_alpha_i) <-
                                  paste0(specification_alpha_i$relation,
                                         "|",
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

## \code{$expand_specification_psi()} expand the specification table for covariances. ##
lslxModel$set("private",
              "expand_specification_psi",
              function() {
                specification_psi <-
                  do.call(what = rbind.data.frame,
                          args = lapply(
                            X = self$name_group,
                            FUN = function(name_group_i) {
                              if (length(self$name_exogenous) > 1) {
                                relation_psi_i <-
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
                                relation_psi_i <-
                                  setdiff(
                                    x = paste0(self$name_eta,
                                               "<->",
                                               self$name_eta),
                                    y = self$specification$relation[self$specification$group == name_group_i]
                                  )
                              }
                              if (length(relation_psi_i) > 1) {
                                left_psi_i <-
                                  substr(relation_psi_i,
                                         start = 1,
                                         stop = regexpr("<->", relation_psi_i) - 1)
                                right_psi_i <-
                                  substr(
                                    relation_psi_i,
                                    start = regexpr("<->", relation_psi_i) + 3,
                                    stop = nchar(relation_psi_i)
                                  )
                                specification_psi_i <-
                                  data.frame(
                                    relation = relation_psi_i,
                                    left = left_psi_i,
                                    right = right_psi_i,
                                    group = name_group_i,
                                    reference = ifelse(
                                      is.na(self$reference_group),
                                      FALSE,
                                      ifelse(name_group_i == self$reference_group,
                                             TRUE,
                                             FALSE)
                                    ),
                                    matrice = "psi",
                                    block =
                                      ifelse(
                                        left_psi_i %in% self$name_factor,
                                        ifelse(right_psi_i %in% self$name_factor,
                                               "f<->f",
                                               "f<->y"),
                                        ifelse(right_psi_i %in% self$name_factor,
                                               "y<->f",
                                               "y<->y")
                                      ),
                                    type = "free",
                                    start = NA_real_,
                                    stringsAsFactors = FALSE
                                  )
                                rownames(specification_psi_i) <-
                                  paste(specification_psi_i$relation,
                                        specification_psi_i$group,
                                        sep = "|")
                                return(specification_psi_i)
                              }
                            }
                          ))
                self$specification <-
                  rbind(self$specification,
                        specification_psi,
                        stringsAsFactors = FALSE)
              })