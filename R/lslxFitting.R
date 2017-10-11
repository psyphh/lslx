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
                  self$control$lambda_grid <-
                    sort(self$control$lambda_grid, decreasing = TRUE)
                  self$control$gamma_grid <-
                    sort(self$control$gamma_grid, decreasing = TRUE)
                  if (length(data$response) > 0) {
                    self$control$raw <- TRUE
                  } else {
                    self$control$raw <- FALSE
                  }
                  if (self$control$missing_method == "default") {
                    if (self$control$raw) {
                      self$control$missing_method <- "two_stage"
                    } else {
                      self$control$missing_method <- "listwise_deletion"
                    }
                  }
                  if (self$control$start_method == "default") {
                    self$control$start_method <- "MH"
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
                      
                      theta_start = model$specification$start,
                      identity_y = matrix(),
                      identity_eta = matrix(),
                      identity_y2 = Matrix::Matrix(),
                      duplication_y = Matrix::Matrix(),
                      elimination_y = Matrix::Matrix(),
                      duplication_eta = Matrix::Matrix(),
                      commutation_y = Matrix::Matrix()
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
                  self$reduced_model$identity_y <-
                    diag(1, self$reduced_model$n_response)
                  self$reduced_model$identity_eta <-
                    diag(1, self$reduced_model$n_eta)
                  self$reduced_model$identity_y2 <-
                    as(
                      diag(
                        1,
                        self$reduced_model$n_response * self$reduced_model$n_response
                      ),
                      "dgCMatrix"
                    )
                  self$reduced_model$duplication_y <-
                    private$create_duplication_matrice(i = self$reduced_model$n_response)
                  self$reduced_model$elimination_y <-
                    as(
                      solve(
                        Matrix::t(self$reduced_model$duplication_y) %*%
                          self$reduced_model$duplication_y
                      ) %*%
                        Matrix::t(self$reduced_model$duplication_y),
                      "dgCMatrix"
                    )
                  self$reduced_model$duplication_eta <-
                    private$create_duplication_matrice(i = self$reduced_model$n_eta)
                  self$reduced_model$commutation_y <-
                    private$create_commutation_matrice(i = self$reduced_model$n_response)
                })

lslxFitting$set("private",
                "create_commutation_matrice",
                function(i) {
                  commutation_i <- Matrix::Matrix(0, i * i, i * i)
                  idx_1 <- 1:(i * i)
                  idx_2 <- c(t(matrix(c(1:(
                    i * i
                  )), i, i)))
                  idx_3 <- idx_2 - 1
                  idx <- (i * i) * idx_3 + idx_1
                  commutation_i[idx] <- 1
                  return(commutation_i)
                })



lslxFitting$set("private",
                "create_duplication_matrice",
                function(i) {
                  duplication_i <- diag(i)
                  idx <- seq(i * (i + 1) / 2)
                  duplication_i[lower.tri(duplication_i,
                                          diag = TRUE)] <-
                    idx
                  duplication_i[upper.tri(duplication_i)] <-
                    t(duplication_i)[upper.tri(duplication_i)]
                  duplication_i <-
                    outer(c(duplication_i),
                          idx,
                          function(x, y)
                            ifelse(x == y, 1, 0))
                  return(Matrix::Matrix(duplication_i))
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
                  if (self$control$raw) {
                    idc_use <-
                      lapply(
                        X = data$pattern,
                        FUN = function(pattern_i) {
                          idc_use_i <-
                            apply(X = pattern_i,
                                  MARGIN = 1,
                                  FUN = prod)
                          return(as.logical(idc_use_i))
                        }
                      )
                    if (self$control$missing_method == "two_stage") {
                      response <- data$response
                      pattern <- data$pattern
                      weight <- data$weight
                    } else if (self$control$missing_method == "listwise_deletion") {
                      response <- 
                        mapply(FUN = function(response_i,
                                              idc_use_i) {
                          response_i <- 
                            response_i[idc_use_i, , drop = FALSE]
                          return(response_i)
                          },
                          data$response,
                          idc_use,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE)
                      pattern <- 
                        mapply(FUN = function(pattern_i,
                                              idc_use_i) {
                          pattern_i <- 
                            pattern_i[idc_use_i, , drop = FALSE]
                          return(pattern_i)
                        },
                        data$pattern,
                        idc_use,
                        SIMPLIFY = FALSE,
                        USE.NAMES = TRUE)
                      weight <- 
                        mapply(FUN = function(weight_i,
                                              idc_use_i) {
                          weight_i <- weight_i[idc_use_i]
                          weight_i <- weight_i / sum(weight_i)
                          return(weight_i)
                        },
                        data$weight,
                        idc_use,
                        SIMPLIFY = FALSE,
                        USE.NAMES = TRUE)
                    } else {
                    }
                    self$reduced_data$n_observation <-
                      sum(sapply(X = response, FUN = nrow))
                    self$reduced_data$n_complete_observation <- sum(unlist(idc_use))
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
                    self$reduced_data$saturated_cov <-
                      lapply(
                        X = self$reduced_data$saturated_cov,
                        FUN = function(saturated_cov_i) {
                          diag(saturated_cov_i) <-
                            diag(saturated_cov_i) + self$control$ridge_cov
                          return(saturated_cov_i)
                        }
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
                    compute_saturated_moment_acov_cpp(
                      y_obs = y_obs,
                      w = w,
                      m_idx = m_idx,
                      m2_idx = m2_idx,
                      saturated_mean = self$reduced_data$saturated_mean,
                      saturated_cov = self$reduced_data$saturated_cov,
                      saturated_moment_acov = self$reduced_data$saturated_moment_acov
                    )
                    self$reduced_data$n_missing_pattern <- nlevels(unlist(m_factor)) 
                  } else {
                    if (self$control$missing_method == "two_stage") {
                      stop("Argument 'missing_method' cannot be 'two_stage' when only moment data is available.")
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
                      mapply(
                        FUN = function(saturated_cov_i,
                                       sample_proportion_i) {
                          saturated_moment_acov_i <-
                            matrix(0,
                                   self$reduced_model$n_moment,
                                   self$reduced_model$n_moment)
                          saturated_cov_i_inv <-
                            solve(saturated_cov_i)
                          saturated_moment_acov_i[1:self$reduced_model$n_response,
                                                  1:self$reduced_model$n_response] <-
                            saturated_cov_i_inv
                          saturated_moment_acov_i[
                            (self$reduced_model$n_response + 1):self$reduced_model$n_moment,
                            (self$reduced_model$n_response + 1):self$reduced_model$n_moment] <-
                            as.matrix(
                              0.5 * Matrix::t(self$reduced_model$duplication_y) %*%
                                (saturated_cov_i_inv %x% saturated_cov_i_inv) %*%
                                self$reduced_model$duplication_y
                            )
                          saturated_moment_acov_i <-
                            solve(saturated_moment_acov_i) / 
                            (sample_proportion_i * self$reduced_data$n_observation)
                          return(saturated_moment_acov_i)
                        },
                        self$reduced_data$saturated_cov,
                        self$reduced_data$sample_proportion,
                        SIMPLIFY = FALSE,
                        USE.NAMES = TRUE
                      )
                    self$reduced_data$n_missing_pattern <- 1
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
                  if (self$control$start_method == "MH") {
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
                            coefficient_matrice_start[[i]][
                              cbind(self$reduced_model$theta_left_idx[idc_i0],
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
                        length(self$control$gamma_grid)
                    )
                  self$fitted_result$information_criterion <-
                    vector(
                      mode = "list",
                      length = length(self$control$lambda_grid) *
                        length(self$control$gamma_grid)
                    )
                  self$fitted_result$fit_indice <-
                    vector(
                      mode = "list",
                      length = length(self$control$lambda_grid) *
                        length(self$control$gamma_grid)
                    )
                  self$fitted_result$coefficient <-
                    vector(
                      mode = "list",
                      length = length(self$control$lambda_grid) *
                        length(self$control$gamma_grid)
                    )
                })
