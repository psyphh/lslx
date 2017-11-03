lslxFitting <-
  R6::R6Class(
    classname = "lslxFitting",
    public = list(
      control = "list",
      reduced_model = "list",
      reduced_data = "list",
      supplied_result = "list",
      fitted_result = "list"
    )
  )

lslxFitting$set("public",
                "initialize",
                function(model,
                         data,
                         control) {
                  private$initialize_control(model = model,
                                             data = data,
                                             control = control)
                  private$initialize_reduced_model(model = model)
                  private$initialize_reduced_data(data = data)
                  private$initialize_supplied_result()
                  private$initialize_fitted_result()
                })


lslxFitting$set("private",
                "initialize_control",
                function(model,
                         data,
                         control) {
                  self$control <- control
                  if (length(data$response) > 0) {
                    self$control$response <- TRUE
                  } else {
                    self$control$response <- FALSE
                  }
                  if (length(data$auxiliary) > 0) {
                    self$control$auxiliary <- TRUE
                  } else {
                    self$control$auxiliary <- FALSE
                  }
                  if (self$control$penalty_method %in% c("mcp", "lasso")) {
                    self$control$regularizer <- TRUE
                  } else {
                    self$control$regularizer <- FALSE
                  }
                  if (self$control$lambda_grid[[1]] == "default") {
                    self$control$lambda_grid <- 0
                  } else {
                    if (self$control$penalty_method == "none") {
                      self$control$lambda_grid <- 0
                    } else if (self$control$penalty_method == "lasso") {
                      if (any(self$control$lambda_grid < 0)) {
                        stop(
                          "When argument 'penalty_method' is set as 'lasso', any element in argument 'lambda_grid' cannot be smaller than 0."
                        )
                      }
                    } else if (self$control$penalty_method == "mcp") {
                      if (any(self$control$lambda_grid < 0)) {
                        stop(
                          "When argument 'penalty_method' is set as 'mcp', any element in argument 'lambda_grid' cannot be smaller than 0."
                        )
                      }
                    } else {
                      
                    }
                  }
                  if (self$control$delta_grid[[1]] == "default") {
                    self$control$delta_grid <- Inf
                  } else {
                    if (self$control$penalty_method %in% c("none", "lasso")) {
                      self$control$delta_grid <- Inf
                    } else if (self$control$penalty_method == "mcp") {
                      if (any(self$control$delta_grid <= 0)) {
                        stop(
                          "When argument 'penalty_method' is set as 'mcp', any element in argument 'delta_grid' must be positive."
                        )
                      }
                    } else {
                      
                    }
                  }
                  self$control$lambda_grid <-
                    sort(self$control$lambda_grid, decreasing = TRUE)
                  self$control$delta_grid <-
                    sort(self$control$delta_grid, decreasing = TRUE)
                  if (self$control$algorithm == "default") {
                    if (self$control$regularizer) {
                      self$control$algorithm <- "fisher"
                    } else {
                      self$control$algorithm <- "bfgs"
                    }
                  }
                  if (self$control$missing_method == "default") {
                    if (self$control$response) {
                      self$control$missing_method <- "two_stage"
                    } else {
                      self$control$missing_method <- "listwise_deletion"
                    }
                  }
                  if (self$control$start_method == "default") {
                    self$control$start_method <- "mh"
                  }
                })


lslxFitting$set("private",
                "initialize_reduced_model",
                function(model) {
                  self$reduced_model <-
                    list(
                      n_response = length(model$name_response),
                      n_factor =  length(model$name_factor),
                      n_eta = length(model$name_eta),
                      n_moment = length(model$name_response) *
                        (length(model$name_response) + 3) / 2,
                      n_group = length(model$name_group),
                      n_theta = nrow(model$specification),
                      theta_name = rownames(model$specification),
                      theta_matrice_idx =
                        ifelse(
                          model$specification$matrice == "alpha",
                          1L,
                          ifelse(
                            model$specification$matrice == "beta",
                            2L,
                            ifelse(model$specification$matrice == "psi",
                                   3L,
                                   0L)
                          )
                        ),
                      theta_left_idx = match(model$specification$left, model$name_eta),
                      theta_right_idx = ifelse(
                        model$specification$matrice == "alpha",
                        1L,
                        ifelse(
                          model$specification$matrice == "tau",
                          NA_integer_,
                          match(model$specification$right, model$name_eta)
                        )
                      ),
                      theta_flat_idx = NA_integer_,
                      theta_group_idx = ifelse(
                        rep(
                          is.na(model$reference_group),
                          length(model$specification$group)
                        ),
                        match(model$specification$group,
                              model$name_group),
                        ifelse(
                          model$specification$group ==
                            model$reference_group,
                          0L,
                          match(model$specification$group,
                                model$name_group)
                        )
                      ),
                      theta_is_free = (model$specification$type == "free"),
                      theta_is_pen = (model$specification$type == "pen"),
                      theta_is_diag =
                        ifelse(
                          model$specification$matrice == "psi" &
                            (model$specification$left ==
                               model$specification$right),
                          TRUE,
                          FALSE
                        ),
                      theta_start = model$specification$start
                    )
                  self$reduced_model$theta_flat_idx <-
                    ifelse(
                      self$reduced_model$theta_matrice_idx == 1,
                      self$reduced_model$theta_left_idx,
                      ifelse(
                        self$reduced_model$theta_matrice_idx == 2,
                        self$reduced_model$n_eta *
                          (self$reduced_model$theta_right_idx - 1L) +
                          self$reduced_model$theta_left_idx,
                        ifelse(
                          self$reduced_model$theta_matrice_idx == 3,
                          as.integer(
                            self$reduced_model$n_eta *
                              (self$reduced_model$theta_right_idx - 1L) +
                              self$reduced_model$theta_left_idx -
                              self$reduced_model$theta_right_idx *
                              (self$reduced_model$theta_right_idx - 1L) / 2L
                          ),
                          NA_integer_
                        )
                      )
                    )
                })


lslxFitting$set("private",
                "initialize_reduced_data",
                function(data) {
                  self$reduced_data <-
                    list(
                      n_observation = integer(),
                      n_complete_observation = integer(),
                      n_missing_pattern = integer(),
                      sample_proportion = list(),
                      saturated_cov = list(),
                      saturated_mean = list(),
                      saturated_moment_acov = list()
                    )
                  if (self$control$response) {
                    idc_complete <-
                      lapply(
                        X = data$pattern,
                        FUN = function(pattern_i) {
                          idc_complete_i <-
                            apply(X = pattern_i,
                                  MARGIN = 1,
                                  FUN = prod)
                          return(as.logical(idc_complete_i))
                        }
                      )
                    if (self$control$missing_method == "two_stage") {
                      idc_use <-
                        lapply(
                          X = data$pattern,
                          FUN = function(pattern_i) {
                            return(rowSums(pattern_i) > 0)
                          }
                        )
                      if (self$control$auxiliary) {
                        response <-
                          mapply(
                            FUN = function(response_i,
                                           auxiliary_i,
                                           idc_use_i) {
                              return(cbind(response_i[idc_use_i, , drop = FALSE],
                                           auxiliary_i[idc_use_i, , drop = FALSE]))
                            },
                            data$response,
                            data$auxiliary,
                            idc_use,
                            SIMPLIFY = FALSE,
                            USE.NAMES = TRUE
                          )
                        pattern <-
                          mapply(
                            FUN = function(pattern_i,
                                           auxiliary_i,
                                           idc_use_i) {
                              return(cbind(pattern_i[idc_use_i, , drop = FALSE],!is.na(auxiliary_i[idc_use_i, , drop = FALSE])))
                            },
                            data$pattern,
                            data$auxiliary,
                            idc_use,
                            SIMPLIFY = FALSE,
                            USE.NAMES = TRUE
                          )
                      } else {
                        response <-
                          mapply(
                            FUN = function(response_i,
                                           idc_use_i) {
                              return(response_i[idc_use_i, , drop = FALSE])
                            },
                            data$response,
                            idc_use,
                            SIMPLIFY = FALSE,
                            USE.NAMES = TRUE
                          )
                        pattern <-
                          mapply(
                            FUN = function(pattern_i,
                                           idc_use_i) {
                              return(pattern_i[idc_use_i, , drop = FALSE])
                            },
                            data$pattern,
                            idc_use,
                            SIMPLIFY = FALSE,
                            USE.NAMES = TRUE
                          )
                      }
                      weight <-
                        mapply(
                          FUN = function(weight_i,
                                         idc_use_i) {
                            weight_i <- weight_i[idc_use_i, ]
                            weight_i <- weight_i / sum(weight_i)
                            return(weight_i)
                          },
                          data$weight,
                          idc_use,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                    } else if (self$control$missing_method == "listwise_deletion") {
                      response <-
                        mapply(
                          FUN = function(response_i,
                                         idc_complete_i) {
                            return(response_i[idc_complete_i, , drop = FALSE])
                          },
                          data$response,
                          idc_complete,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                      pattern <-
                        mapply(
                          FUN = function(pattern_i,
                                         idc_complete_i) {
                            return(pattern_i[idc_complete_i, , drop = FALSE])
                          },
                          data$pattern,
                          idc_complete,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                      weight <-
                        mapply(
                          FUN = function(weight_i,
                                         idc_complete_i) {
                            weight_i <- weight_i[idc_complete_i, ]
                            weight_i <- weight_i / sum(weight_i)
                            return(weight_i)
                          },
                          data$weight,
                          idc_complete,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                    } else {
                      
                    }
                    self$reduced_data$n_observation <-
                      sum(sapply(X = response, FUN = nrow))
                    self$reduced_data$n_complete_observation <-
                      sum(unlist(idc_complete))
                    self$reduced_data$sample_proportion <-
                      lapply(
                        X = lapply(X = response, FUN = nrow),
                        FUN = function(sample_size_i) {
                          sample_proportion_i <-
                            sample_size_i / self$reduced_data$n_observation
                          return(sample_proportion_i)
                        }
                      )
                    self$reduced_data$saturated_mean <-
                      lapply(
                        X = response,
                        FUN = function(response_i) {
                          as.matrix(colMeans(response_i, na.rm = TRUE))
                        }
                      )
                    self$reduced_data$saturated_cov <-
                      lapply(X = response,
                             FUN = cov,
                             use = "pairwise.complete.obs")
                    m_factor <-
                      lapply(
                        X = pattern,
                        FUN = function(pattern_i) {
                          as.factor(apply(
                            X = pattern_i,
                            MARGIN = 1,
                            FUN = function(pattern_ij) {
                              do.call(what = paste,
                                      args = as.list(c(pattern_ij, sep = "/")))
                            }
                          ))
                        }
                      )
                    m_idx <-
                      lapply(
                        X = m_factor,
                        FUN = function(m_factor_i) {
                          m_idx_i <-
                            lapply(
                              X = strsplit(levels(m_factor_i),
                                           split = "/"),
                              FUN = function(x) {
                                return(as.integer(which(as.logical(x))) - 1)
                              }
                            )
                          names(m_idx_i) <- levels(m_factor_i)
                          return(m_idx_i)
                        }
                      )
                    y_obs <-
                      mapply(
                        FUN = function(response_i,
                                       m_factor_i,
                                       m_idx_i) {
                          y_i <- split(response_i, m_factor_i)
                          mapply(
                            FUN = function(y_ij,
                                           m_idx_ij) {
                              y_obs_ij <- as.matrix(y_ij[, (m_idx_ij + 1) , drop = FALSE])
                              storage.mode(y_obs_ij) <- "double"
                              return(y_obs_ij)
                            },
                            y_i,
                            m_idx_i,
                            SIMPLIFY = FALSE,
                            USE.NAMES = TRUE
                          )
                        },
                        response,
                        m_factor,
                        m_idx,
                        SIMPLIFY = FALSE,
                        USE.NAMES = TRUE
                      )
                    w <-
                      mapply(
                        FUN = function(weight_i,
                                       m_factor_i) {
                          split(weight_i, m_factor_i)
                        },
                        weight,
                        m_factor,
                        SIMPLIFY = FALSE,
                        USE.NAMES = TRUE
                      )
                    
                    compute_saturated_moment_cpp(
                      y_obs = y_obs,
                      w = w,
                      m_idx = m_idx,
                      saturated_mean = self$reduced_data$saturated_mean,
                      saturated_cov = self$reduced_data$saturated_cov,
                      iter_other_max = self$control$iter_other_max,
                      tol_other = self$control$tol_other
                    )
                    m2_idx <-
                      lapply(
                        X = m_factor,
                        FUN = function(m_factor_i) {
                          m2_idx_i <-
                            lapply(
                              X = strsplit(levels(m_factor_i),
                                           split = "/"),
                              FUN = function(x) {
                                m2_idx_i <- tcrossprod(as.matrix(as.logical(x)))
                                m2_idx_i <-
                                  as.logical(m2_idx_i[lower.tri(m2_idx_i, diag = TRUE)])
                                return(as.integer(which(m2_idx_i)) - 1)
                              }
                            )
                          names(m2_idx_i) <- levels(m_factor_i)
                          return(m2_idx_i)
                        }
                      )
                    self$reduced_data$saturated_moment_acov <-
                      lapply(
                        X = self$reduced_data$saturated_cov,
                        FUN = function(saturated_cov_i) {
                          
                        }
                      )
                    compute_saturated_moment_acov_response_cpp(
                      y_obs = y_obs,
                      w = w,
                      m_idx = m_idx,
                      m2_idx = m2_idx,
                      saturated_mean = self$reduced_data$saturated_mean,
                      saturated_cov = self$reduced_data$saturated_cov,
                      saturated_moment_acov = self$reduced_data$saturated_moment_acov
                    )
                    self$reduced_data$n_missing_pattern <-
                      nlevels(unlist(m_factor))
                  } else {
                    if (self$control$missing_method == "two_stage") {
                      stop(
                        "Argument 'missing_method' cannot be 'two_stage' when only moment data is available."
                      )
                    }
                    self$reduced_data$n_observation <-
                      sum(unlist(data$sample_size))
                    self$reduced_data$n_complete_observation <-
                      self$reduced_data$n_observation
                    self$reduced_data$sample_proportion <-
                      lapply(
                        X = data$sample_size,
                        FUN = function(sample_size_i) {
                          sample_proportion_i <-
                            sample_size_i / self$reduced_data$n_observation
                          return(sample_proportion_i)
                        }
                      )
                    self$reduced_data$saturated_cov <-
                      lapply(
                        X = data$sample_cov,
                        FUN = function(sample_cov_i) {
                          diag(sample_cov_i) <-
                            diag(sample_cov_i) + self$control$ridge_cov
                          return(sample_cov_i)
                        }
                      )
                    self$reduced_data$saturated_mean <-
                      lapply(
                        X = data$sample_mean,
                        FUN = function(sample_mean_i) {
                          sample_mean_i <- as.matrix(sample_mean_i)
                          return(sample_mean_i)
                        }
                      )
                    self$reduced_data$saturated_moment_acov <-
                      lapply(
                        X = self$reduced_data$saturated_cov,
                        FUN = function(saturated_cov_i) {
                          
                        }
                      )
                    compute_saturated_moment_acov_moment_cpp(
                      n_observation = self$reduced_data$n_observation,
                      sample_proportion = self$reduced_data$sample_proportion,
                      saturated_cov = self$reduced_data$saturated_cov,
                      saturated_moment_acov = self$reduced_data$saturated_moment_acov
                    )
                    self$reduced_data$n_missing_pattern <- 1
                  }
                  
                  if (self$control$auxiliary) {
                    y_name <-
                      c(colnames(data$response[[1]]), colnames(data$auxiliary[[1]]))
                  } else {
                    if (self$control$response) {
                      y_name <- colnames(data$response[[1]])
                    } else {
                      y_name <- colnames(data$sample_cov[[1]])
                    }
                  }
                  y2_name <- outer(
                    y_name,
                    y_name,
                    FUN = function(y_name_i, y_name_j) {
                      return(paste(y_name_i, y_name_j, sep = "*"))
                    }
                  )
                  y2_name <-
                    y2_name[lower.tri(y2_name, diag = TRUE)]
                  self$reduced_data$saturated_mean <-
                    lapply(
                      X = self$reduced_data$saturated_mean,
                      FUN = function(saturated_mean_i) {
                        rownames(saturated_mean_i) <- y_name
                        return(saturated_mean_i)
                      }
                    )
                  self$reduced_data$saturated_cov <-
                    lapply(
                      X = self$reduced_data$saturated_cov,
                      FUN = function(saturated_cov_i) {
                        diag(saturated_cov_i) <-
                          diag(saturated_cov_i) + self$control$ridge_cov
                        rownames(saturated_cov_i) <- y_name
                        colnames(saturated_cov_i) <- y_name
                        return(saturated_cov_i)
                      }
                    )
                  self$reduced_data$saturated_moment_acov <-
                    lapply(
                      X = self$reduced_data$saturated_moment_acov,
                      FUN = function(saturated_moment_acov_i) {
                        rownames(saturated_moment_acov_i) <-
                          c(y_name, y2_name)
                        colnames(saturated_moment_acov_i) <-
                          c(y_name, y2_name)
                        return(saturated_moment_acov_i)
                      }
                    )
                  if (self$control$auxiliary) {
                    y_name <- colnames(data$response[[1]])
                    y2_name <- outer(
                      y_name,
                      y_name,
                      FUN = function(y_name_i, y_name_j) {
                        return(paste(y_name_i, y_name_j, sep = "*"))
                      }
                    )
                    y2_name <-
                      y2_name[lower.tri(y2_name, diag = TRUE)]
                    self$reduced_data$saturated_mean <-
                      lapply(
                        X = self$reduced_data$saturated_mean,
                        FUN = function(saturated_mean_i) {
                          return(saturated_mean_i[y_name, , drop = FALSE])
                        }
                      )
                    self$reduced_data$saturated_cov <-
                      lapply(
                        X = self$reduced_data$saturated_cov,
                        FUN = function(saturated_cov_i) {
                          return(saturated_cov_i[y_name, y_name, drop = FALSE])
                        }
                      )
                    self$reduced_data$saturated_moment_acov <-
                      lapply(
                        X = self$reduced_data$saturated_moment_acov,
                        FUN = function(saturated_moment_acov_i) {
                          return(saturated_moment_acov_i[c(y_name, y2_name),
                                                         c(y_name, y2_name),
                                                         drop = FALSE])
                        }
                      )
                  } else {
                    
                  }
                })


lslxFitting$set("private",
                "initialize_supplied_result",
                function() {
                  self$supplied_result <- list()
                  private$compute_baseline_model()
                  private$compute_saturated_model()
                  private$compute_fitted_start()
                })


lslxFitting$set("private",
                "compute_fitted_start",
                function() {
                  if (self$control$start_method == "mh") {
                    saturated_cov_pool <-
                      Reduce(
                        "+",
                        mapply(
                          FUN = function(sample_proportion_i,
                                         saturated_cov_i) {
                            return(sample_proportion_i * saturated_cov_i)
                          },
                          self$reduced_data$sample_proportion,
                          self$reduced_data$saturated_cov,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                      )
                    saturated_mean_pool <-
                      Reduce(
                        "+",
                        mapply(
                          FUN = function(sample_proportion_i,
                                         saturated_mean_i) {
                            return(sample_proportion_i * saturated_mean_i)
                          },
                          self$reduced_data$sample_proportion,
                          self$reduced_data$saturated_mean,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                      )
                    idc_beta <-
                      (self$reduced_model$theta_matrice_idx == 2 &
                         self$reduced_model$theta_is_free) |
                      (
                        self$reduced_model$theta_matrice_idx == 2 &
                          (
                            !self$reduced_model$theta_is_free &
                              !self$reduced_model$theta_is_pen
                          ) &
                          (self$reduced_model$theta_start != 0)
                      ) |
                      (
                        self$reduced_model$theta_matrice_idx == 2 &
                          self$reduced_model$theta_is_pen &
                          !(
                            self$reduced_model$theta_left_idx <=
                              self$reduced_model$n_response &
                              self$reduced_model$theta_right_idx >
                              self$reduced_model$n_response
                          )
                      )
                    idx_beta <-
                      strsplit(x = unique(
                        paste0(
                          self$reduced_model$theta_left_idx[idc_beta],
                          ",",
                          self$reduced_model$theta_right_idx[idc_beta]
                        )
                      ),
                      split = ",")
                    theta_left_idx_beta <-
                      sapply(
                        X = idx_beta,
                        FUN = function(idx_beta_i) {
                          as.integer(idx_beta_i[1])
                        }
                      )
                    theta_right_idx_beta <-
                      sapply(
                        X = idx_beta,
                        FUN = function(idx_beta_i) {
                          as.integer(idx_beta_i[2])
                        }
                      )
                    if (self$reduced_model$n_factor > 0) {
                      cov_eta <-
                        matrix(0,
                               self$reduced_model$n_eta,
                               self$reduced_model$n_eta)
                      cov_eta[1:self$reduced_model$n_response,
                              1:self$reduced_model$n_response] <-
                        cov2cor(saturated_cov_pool)
                      for (i in (self$reduced_model$n_response + 1):(self$reduced_model$n_eta)) {
                        theta_left_idx_beta_i <-
                          theta_left_idx_beta[theta_right_idx_beta == i &
                                                theta_left_idx_beta <= self$reduced_model$n_response]
                        cov_sum_i <-
                          sum(cov_eta[theta_left_idx_beta_i,
                                      theta_left_idx_beta_i])
                        for (j in seq_len(i)) {
                          if (j <= self$reduced_model$n_response) {
                            cov_eta_ij <-
                              sum(cov_eta[j, theta_left_idx_beta_i]) / sqrt(abs(cov_sum_i))
                            cov_eta[i, j] <- cov_eta_ij
                            cov_eta[j, i] <- cov_eta_ij
                          } else if (j < i) {
                            theta_left_idx_beta_j <-
                              theta_left_idx_beta[theta_right_idx_beta == j &
                                                    theta_left_idx_beta <= self$reduced_model$n_response]
                            cov_sum_j <-
                              sum(cov_eta[theta_left_idx_beta_j,
                                          theta_left_idx_beta_j])
                            cov_eta_ij <-
                              sum(cov_eta[theta_left_idx_beta_j, theta_left_idx_beta_i]) /
                              sqrt(abs(cov_sum_i) * abs(cov_sum_j))
                            cov_eta[i, j] <- cov_eta_ij
                            cov_eta[j, i] <- cov_eta_ij
                          } else {
                            cov_eta[j, i] <- 1
                          }
                        }
                      }
                      cov_eta <-
                        diag(c(sqrt(diag(
                          saturated_cov_pool
                        )),
                        rep(1, self$reduced_model$n_factor))) %*%
                        cov_eta %*%
                        diag(c(sqrt(diag(
                          saturated_cov_pool
                        )),
                        rep(1, self$reduced_model$n_factor)))
                    } else {
                      cov_eta <- saturated_cov_pool
                    }
                    beta_start <-
                      matrix(0,
                             self$reduced_model$n_eta,
                             self$reduced_model$n_eta)
                    for (i in seq_len(self$reduced_model$n_eta)) {
                      theta_right_idx_beta_i <-
                        theta_right_idx_beta[theta_left_idx_beta == i]
                      if (length(theta_right_idx_beta_i) > 0) {
                        beta_start_i <-
                          solve(cov_eta[theta_right_idx_beta_i,
                                        theta_right_idx_beta_i,
                                        drop = FALSE]) %*%
                          cov_eta[theta_right_idx_beta_i, i, drop = FALSE]
                        beta_start[i, theta_right_idx_beta_i] <-
                          beta_start_i
                      }
                    }
                    idc_psi <-
                      (self$reduced_model$theta_matrice_idx == 3 &
                         self$reduced_model$theta_is_free) |
                      (
                        self$reduced_model$theta_matrice_idx == 3 &
                          (
                            !self$reduced_model$theta_is_free &
                              !self$reduced_model$theta_is_pen
                          ) &
                          (self$reduced_model$theta_start != 0)
                      )
                    idx_psi <-
                      strsplit(x = unique(
                        paste0(
                          self$reduced_model$theta_left_idx[idc_psi],
                          ",",
                          self$reduced_model$theta_right_idx[idc_psi]
                        )
                      ),
                      split = ",")
                    theta_left_idx_psi <-
                      sapply(
                        X = idx_psi,
                        FUN = function(idx_psi_i) {
                          as.integer(idx_psi_i[1])
                        }
                      )
                    theta_right_idx_psi <-
                      sapply(
                        X = idx_psi,
                        FUN = function(idx_psi_i) {
                          as.integer(idx_psi_i[2])
                        }
                      )
                    psi_start <-
                      matrix(0,
                             self$reduced_model$n_eta,
                             self$reduced_model$n_eta)
                    psi_saturated <-
                      (diag(1, self$reduced_model$n_eta) - beta_start) %*%
                      cov_eta %*% t((diag(1, self$reduced_model$n_eta) - beta_start))
                    psi_start[cbind(theta_left_idx_psi, theta_right_idx_psi)] <-
                      psi_saturated[cbind(theta_left_idx_psi, theta_right_idx_psi)]
                    psi_start[cbind(theta_right_idx_psi, theta_left_idx_psi)] <-
                      psi_saturated[cbind(theta_right_idx_psi, theta_left_idx_psi)]
                    idc_alpha <-
                      (self$reduced_model$theta_matrice_idx == 1 &
                         self$reduced_model$theta_is_free) |
                      (
                        self$reduced_model$theta_matrice_idx == 1 &
                          (
                            !self$reduced_model$theta_is_free &
                              !self$reduced_model$theta_is_pen
                          ) &
                          (self$reduced_model$theta_start != 0)
                      )
                    theta_left_idx_alpha <-
                      unique(self$reduced_model$theta_left_idx[idc_alpha])
                    if (self$reduced_model$n_factor > 0) {
                      mean_eta <-
                        matrix(0, self$reduced_model$n_eta, 1)
                      mean_eta[1:self$reduced_model$n_response, 1] <-
                        saturated_mean_pool
                      for (i in (self$reduced_model$n_response + 1):(self$reduced_model$n_eta)) {
                        if (i %in% theta_left_idx_alpha) {
                          theta_left_idx_beta_i <-
                            theta_left_idx_beta[theta_right_idx_beta == i]
                          mean_eta[i, 1] <-
                            solve(t(beta_start[theta_left_idx_beta_i, i, drop = FALSE]) %*%
                                    beta_start[theta_left_idx_beta_i, i, drop = FALSE]) %*%
                            t(beta_start[theta_left_idx_beta_i, i, drop = FALSE]) %*%
                            mean_eta[theta_left_idx_beta_i, 1, drop = FALSE]
                        } else {
                          mean_eta[i, 1] <- 0
                        }
                      }
                    } else {
                      mean_eta <- saturated_mean_pool
                    }
                    alpha_saturated <-
                      (diag(1, self$reduced_model$n_eta) - beta_start) %*% mean_eta
                    alpha_start <-
                      matrix(0, self$reduced_model$n_eta, 1)
                    alpha_start[theta_left_idx_alpha, 1] <-
                      alpha_saturated[theta_left_idx_alpha, 1]
                    
                    coefficient_matrice_start <-
                      list(
                        alpha_start = alpha_start,
                        beta_start = beta_start,
                        psi_start = psi_start
                      )
                    for (i in 1:3) {
                      if (0 %in% unique(self$reduced_model$theta_group_idx)) {
                        idc_i0 <-
                          (self$reduced_model$theta_matrice_idx == i) &
                          (self$reduced_model$theta_group_idx == 0)
                        
                        self$supplied_result$fitted_start[idc_i0] <-
                          ifelse(
                            is.na(self$reduced_model$theta_start[idc_i0]),
                            coefficient_matrice_start[[i]][cbind(self$reduced_model$theta_left_idx[idc_i0],
                                                                 self$reduced_model$theta_right_idx[idc_i0])],
                            self$reduced_model$theta_start[idc_i0]
                          )
                        for (j in setdiff(unique(self$reduced_model$theta_group_idx), 0)) {
                          idc_ij <-
                            (self$reduced_model$theta_matrice_idx == i) &
                            (self$reduced_model$theta_group_idx == j)
                          self$supplied_result$fitted_start[idc_ij] <-
                            ifelse(
                              is.na(self$reduced_model$theta_start[idc_ij]),
                              0,
                              self$reduced_model$theta_start[idc_ij]
                            )
                        }
                      } else {
                        for (j in seq_len(self$reduced_model$n_group)) {
                          idc_ij <-
                            (self$reduced_model$theta_matrice_idx == i) &
                            (self$reduced_model$theta_group_idx == j)
                          self$supplied_result$fitted_start[idc_ij] <-
                            ifelse(
                              is.na(self$reduced_model$theta_start[idc_ij]),
                              coefficient_matrice_start[[i]][cbind(
                                self$reduced_model$theta_left_idx[idc_ij],
                                self$reduced_model$theta_right_idx[idc_ij]
                              )],
                              self$reduced_model$theta_start[idc_ij]
                            )
                        }
                      }
                    }
                  } else if (self$control$start_method == "heuristic") {
                    self$supplied_result$fitted_start <-
                      ifelse(
                        !is.na(self$reduced_model$theta_start),
                        self$reduced_model$theta_start,
                        ifelse(self$reduced_model$matrice == 2,
                               1,
                               ifelse((self$reduced_model$matrice == 3) &
                                        (
                                          self$reduced_model$theta_left_idx ==
                                            self$reduced_model$theta_right_idx
                                        ),
                                      .05,
                                      0
                               ))
                      )
                    if (any(self$reduced_model$theta_group_idx == 0)) {
                      self$supplied_result$fitted_start <-
                        ifelse(
                          self$reduced_model$theta_group_idx == 0,
                          self$supplied_result$fitted_start,
                          0
                        )
                    }
                  } else {
                    
                  }
                  names(self$supplied_result$fitted_start) <-
                    self$reduced_model$theta_name
                })


lslxFitting$set("private",
                "compute_baseline_model",
                function() {
                  self$supplied_result$baseline_model <- c(
                    loss_value =
                      sum(
                        mapply(
                          FUN = function(saturated_mean_i,
                                         saturated_cov_i,
                                         sample_proportion_i) {
                            baseline_loss_value_i <-
                              sample_proportion_i *
                              (sum(diag(saturated_cov_i %*%
                                          diag(
                                            1 / diag(saturated_cov_i)
                                          ))) -
                                 log(det(saturated_cov_i %*% diag(
                                   1 / diag(saturated_cov_i)
                                 ))) -
                                 self$reduced_model$n_response)
                            return(baseline_loss_value_i)
                          },
                          self$reduced_data$saturated_mean,
                          self$reduced_data$saturated_cov,
                          self$reduced_data$sample_proportion,
                          SIMPLIFY = TRUE
                        )
                      ),
                    n_nonzero_coefficient =
                      self$reduced_model$n_group *
                      (2L * self$reduced_model$n_response),
                    degree_of_freedom =
                      self$reduced_model$n_group *
                      (
                        self$reduced_model$n_moment -
                          2L * self$reduced_model$n_response
                      )
                  )
                })


lslxFitting$set("private",
                "compute_saturated_model",
                function() {
                  self$supplied_result$saturated_model <- c(
                    loss_value = 0,
                    n_nonzero_coefficient = self$reduced_model$n_group *
                      self$reduced_model$n_moment,
                    degree_of_freedom = 0
                  )
                })


lslxFitting$set("private",
                "initialize_fitted_result",
                function() {
                  self$fitted_result <- list()
                  self$fitted_result$numerical_condition <-
                    vector(
                      mode = "list",
                      length = length(self$control$lambda_grid) *
                        length(self$control$delta_grid)
                    )
                  self$fitted_result$information_criterion <-
                    vector(
                      mode = "list",
                      length = length(self$control$lambda_grid) *
                        length(self$control$delta_grid)
                    )
                  self$fitted_result$fit_indice <-
                    vector(
                      mode = "list",
                      length = length(self$control$lambda_grid) *
                        length(self$control$delta_grid)
                    )
                  self$fitted_result$coefficient <-
                    vector(
                      mode = "list",
                      length = length(self$control$lambda_grid) *
                        length(self$control$delta_grid)
                    )
                })
