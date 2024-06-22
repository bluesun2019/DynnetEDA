# library(fda)
# # library(rTensor)
# library(doParallel)
# library(doSNOW)
# library(foreach)

# the main functions to perform node embeddings.
#' @title Perform network embeddings.
#' @description Transform the dynamic network into d-dimensional latent functions for each node via Accelerated fASE.
#'
#'
#'
#' @param adjacency_tensor an n*n*T array where each slice is a snapshot of the dynamic network.
#' @param embedding_dim scalar. the dimension of latent functions d.
#' @param spline_no scalar. the number of B-spline basis q.
#' @param spline_order (optional) scalar. the polynomial order of B-spline basis, default 4
#' @param batch_size (optional) scalar. the size of batches in mini-batch algorithm, default max(50,floor(n/20))
#' @param kernel_scale (optional) scalar. the scaled parameter in the kernel used for smoothing, default two times the quartile of the time intervals.
#' @param timestamp_vec (optional) vector. the observed time points of observation, default 1:T.
#' @param scalable_method  (optional) logical. whether to use the random-projection based SVD and mini-batch methods, default TRUE.
#' @param scalable_dim (optional) scalar. the dimension of random projection if scalable_method = TRUE, default 20.
#' @param scalable_power (optional) scalar. the time of power iterations in power method when using  the random-projection based SVD, default 6.
#' @param update_method (optional) string. "Adam" if using the Adam algorithm for update; "GD" if using the classical gradient descent algorithm, default "Adam".
#' @param adam_parameter (optional) three-element list. a list containing beta1, beta2 and epsilon in the Adam algorithm if update_method = "Adam". default beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8.
#' @param iteration_step (optional) scalar. the number of maximal iteration steps for the gradient algorithm, default 2500.
#' @param epsilon (optional) scalar. the stopping criterion for the gradient algorithm, default 1e-6.
#' @param step_size (optional) scalar. the initial step size of the gradient algorithm, default 0.1.
#' @param epoch_length (optional) scalar. iteration numbers in one loop before calculating the objective function, default 1.
#' @param parallel (optional) logical. whether to parallel when tuning hyperparameters, defult TRUE.
#'
#' @return a list of five elements containing: n*d fd format functions, n*d*q coefficient tensor, corresponding functional basis object, values of the basis object on timestamps of snapshots
#'
#' @export
#' @examples
#' ER_generation <- function(n, p) {
#'   adjacency_matrix <- matrix(0, nrow = n, ncol = n)
#'   for (i in 2:n) {
#'     for (j in 1:(i - 1)) {
#'       adjacency_matrix[i, j] <- sample(c(0, 1), size = 1, prob = c(1 - p, p))
#'     }
#'     adjacency_matrix[i, i] <- 0.5
#'   }
#'   adjacency_matrix + t(adjacency_matrix)
#' }
#' set.seed(10)
#' n_nodes <- 10
#' T <- 10
#' p <- 0.2
#' dynamic_adjacency <- c()
#' dynamic_network_adjacency <- array(0, dim = c(n_nodes, n_nodes, T))
#' for (t in 1:T) {
#'   adj <- ER_generation(n_nodes, p)
#'   dynamic_adjacency <- c(dynamic_adjacency, list(list(adj)))
#'   dynamic_network_adjacency[, , t] <- adj
#' }
#' result_ER1 <- fASE(dynamic_network_adjacency, 2, 4, eps = 1, batch_size = NULL)
#'
fASE <- function(adjacency_tensor, embedding_dim, spline_no, spline_order = 4, batch_size = NULL, kernel_scale = FALSE, timestamp_vec = NULL, scalable_method = TRUE,
                 scalable_dim = 20, scalable_power = 6, update_method = "Adam", adam_parameter = NULL, iteration_step = 2500, epsilon = 1e-6, step_size = 0.1, epoch_length = 1, parallel = FALSE) {
  set.seed(1)

  n_nodes <- dim(adjacency_tensor)[1]
  m_layers <- dim(adjacency_tensor)[3]
  adam_parameter_init <- adam_parameter

  if (is.null(timestamp_vec)) {
    timestamp_vec <- 1:m_layers
  }

  if (spline_no - spline_order + 1 > m_layers) {
    cat("the number of knots should be less than the number of layers!")
    break
  }

  if (scalable_method) {
    if (is.null(batch_size)) {
      if (n_nodes < 50) {
        batch_size <- n_nodes
      } else {
        batch_size <- max(50, floor(n_nodes / 20))
      }
    }
  }

  if (is.null(timestamp_vec)) {
    knots <- seq(1, m_layers, by = (m_layers - 1) / (spline_no - spline_order + 1))
  } else {
    knots <- seq(min(timestamp_vec), max(timestamp_vec), by = (max(timestamp_vec) - min(timestamp_vec)) / (spline_no - spline_order + 1))
  }

  # functions for initialization, calculation of objective functions and gradients.
  fx.hat <- function(x, z, h) {
    dnorm((z - x) / h) / h
  }
  smoothing <- function(y, a) {
    l <- length(a)
    new_mat <- matrix(0, nrow = dim(y)[1], ncol = dim(y)[2])
    for (i in 1:l) {
      new_mat <- new_mat + y[, , i] * a[i]
    }
    new_mat
  }
  KSMOOTH <- function(h, y, x, z) {
    m <- length(x)
    s.hat <- array(0, dim = c(dim(y)[1], dim(y)[2], m))
    for (i in 1:m) {
      a <- fx.hat(x[i], z, h)
      s.hat[, , i] <- smoothing(y, a)
    }
    s.hat
  }
  half_decompose <- function(X, d, s, q, scalable_method) {
    if (scalable_method) {
      set.seed(1)

      n <- dim(X)[1]

      Omega <- matrix(rnorm(n * s), nrow = n, ncol = s)
      A <- Omega

      for (i in 1:q) {
        if (i %% 2 == 1) {
          A <- X %*% A
        } else {
          A <- t(X) %*% A
        }
      }

      qr_result <- qr(A)
      rm(A)
      gc()
      Q <- qr.Q(qr_result)

      X_approx <- t(Q) %*% X %*% Q
      X_approx[which(X_approx != t(X_approx))] <- 0
      eigen_result <- eigen(X_approx)
      rm(X, X_approx)
      gc()

      eigenval <- eigen_result$values[1:d]
      eigenval[eigenval < 0] <- 0

      result <- t(t(Q %*% eigen_result$vectors[, 1:d]) * sqrt(eigenval))
    } else {
      eigen_result <- eigen(X)
      rm(X)
      gc()

      eigenval <- eigen_result$values[1:d]
      eigenval[eigenval < 0] <- 0

      result <- t(t(eigen_result$vectors[, 1:d]) * sqrt(eigenval))
    }

    result
  }
  kernel_embedding_initial <- function(adjacency_tensor, n_nodes, m_layers, timestamp_vec, embedding_dim, spline_no, B_splines, scalable_method, random_dim, mul_iter, h, kernel = "KS", spline_order = 4) {
    grid_no <- spline_no - spline_order + 2

    grids <- seq(min(timestamp_vec), max(timestamp_vec), by = (m_layers - 1) / (grid_no - 1))
    B <- t(fda::eval.basis(grids, B_splines))

    if (kernel == "KS") {
      if (!h) {
        h <- 2 * quantile(diff(timestamp_vec), 0.25)
      }
      K <- KSMOOTH(h, adjacency_tensor, grids, timestamp_vec)
    }
    if (is.null(random_dim)) {
      random_dim <- min(max(100, 10 * embedding_dim), n_nodes)
    }
    if (is.null(mul_iter)) {
      mul_iter <- 5
    }

    Z_init <- array(0, dim = c(n_nodes, embedding_dim, grid_no))

    for (l in 1:grid_no) {
      Z_init2 <- half_decompose(K[, , l], embedding_dim, random_dim, mul_iter, scalable_method)

      if (l > 1) {
        svd_result <- svd(t(Z_init2) %*% Z_init[, , l - 1])
        rotation_matrix <- svd_result$v %*% t(svd_result$u)
        Z_init[, , l] <- Z_init2 %*% rotation_matrix
      } else {
        Z_init[, , l] <- Z_init2
      }
    }

    W_init <- apply(Z_init, c(1, 2), crossprod, solve(t(B) %*% B, t(B)))
    W_init <- aperm(W_init, c(2, 3, 1))
    W_init
  }
  gradient_computation <- function(adjacency_tensor, W, B, BB, r, method = "GD", mini_index = NULL, adam_par = NULL) {
    m <- dim(adjacency_tensor)[3]
    n <- dim(adjacency_tensor)[1]
    d <- dim(W)[2]
    if (is.null(mini_index)) {
      mini_index <- 1:n
    }
    b <- length(mini_index)

    sum_matrices <- function(x, y) {
      x + y
    }


    if (method == "Adam") {
      # parameters of Adam are:
      # beta1, beta2, iter, m_(t-1), v_(t-1), epsilon
      beta1 <- adam_par[[1]]
      beta2 <- adam_par[[2]]
      iter <- adam_par[[3]]
      momentum <- adam_par[[4]]
      variance <- adam_par[[5]]
      epsilon <- adam_par[[6]]

      derivate <- 0
      WW <- array(apply(W, 2, crossprod, W[, r, ]), dim = c(dim(W)[3], dim(W)[3], dim(W)[2]))

      if (b <= 200) {
        for (i in 1:m) {
          temp <- adjacency_tensor[, , i]
          derivate <- derivate - n / b * temp[, mini_index] %*% ((W[, r, ] %*% BB[, , i])[mini_index, ])
        }
        for (i in 1:m) {
          for (s in 1:d) {
            derivate <- derivate + tcrossprod(W[, s, ], t(BB[, , i] %*% WW[, , s] %*% BB[, , i]))
          }
        }
      } else {
        for (i in 1:m) {
          temp <- adjacency_tensor[, , i]
          no_of_subbatches <- ceiling(b / 200)
          for (s in 1:no_of_subbatches) {
            mini_subindex <- (200 * (s - 1) + 1):min(200 * s, b)
            derivate <- derivate - n / b * temp[, mini_index[mini_subindex]] %*% ((W[, r, ] %*% BB[, , i])[mini_index[mini_subindex], ])
          }
        }
        for (i in 1:m) {
          for (s in 1:d) {
            derivate <- derivate + tcrossprod(W[, s, ], t(BB[, , i] %*% WW[, , s] %*% BB[, , i]))
          }
        }
      }


      # derivate = foreach(i=1:m, .combine = "sum_matrices")%dopar%{
      #   temp = adjacency_tensor[,,i]
      #   derivate = - n/b * temp[,mini_index]%*%((W[,r,]%*%B[,i]%*%t(B[,i]))[mini_index,])
      #   for(s in 1:d){
      #    derivate = derivate + W[,s,]%*%B[,i]%*%t(B[,i])%*%(t(W[,s,])%*%W[,r,])%*%B[,i]%*%t(B[,i])
      #   }
      #  derivate
      # }

      momentum <- beta1 * momentum + (1 - beta1) * derivate
      variance <- beta2 * variance + (1 - beta2) * derivate^2
      alpha_now <- sqrt(1 - beta2^iter) / (1 - beta1^iter)
      derivate <- alpha_now * momentum / (sqrt(variance) + epsilon)

      result <- list(derivate, momentum, variance)
    } else if (method == "GD") {
      derivate <- 0
      batch_size <- length(mini_index)

      WW <- array(apply(W, 2, crossprod, W[, r, ]), dim = c(dim(W)[3], dim(W)[3], dim(W)[2]))

      for (i in 1:m) {
        temp <- adjacency_tensor[, , i]
        no_of_batches <- ceiling(dim(adjacency_tensor)[1] / batch_size)
        for (s in 1:no_of_batches) {
          mini_index <- (batch_size * (s - 1) + 1):min(batch_size * s, dim(adjacency_tensor)[1])
          derivate <- derivate - temp[, mini_index] %*% ((W[, r, ] %*% B[, i] %*% t(B[, i]))[mini_index, ])
        }
        for (s in 1:d) {
          derivate <- derivate + W[, s, ] %*% B[, i] %*% t(B[, i]) %*% WW[, , s] %*% B[, i] %*% t(B[, i])
        }
      }

      result <- derivate
    }

    result
  }
  objective_computation <- function(adjacency_tensor, W, B) {
    objective_result <- 0

    m <- dim(adjacency_tensor)[3]
    n <- dim(adjacency_tensor)[1]
    d <- dim(W)[2]

    for (i in 1:m) {
      temp <- adjacency_tensor[, , i]
      for (s in 1:d) {
        tempp <- W[, s, ] %*% B[, i]
        temp <- temp - tempp %*% t(tempp)
        rm(tempp)
        gc()
      }

      objective_result <- objective_result + norm(temp, type = "F")^2
      rm(temp)
      gc()
    }

    objective_result
  }


  B_splines <- fda::create.bspline.basis(rangeval = c(knots[1], knots[length(knots)]), nbasis = spline_no, norder = spline_order, breaks = knots)
  B <- t(fda::eval.basis(timestamp_vec, B_splines))
  BB <- array(apply(B, 2, tcrossprod), dim = c(dim(B)[1], dim(B)[1], dim(B)[2]))

  W_init <- kernel_embedding_initial(adjacency_tensor, n_nodes, m_layers, timestamp_vec, embedding_dim, spline_no, B_splines, scalable_method, scalable_dim, scalable_power, kernel_scale, spline_order = spline_order)
  W <- W_init

  if (scalable_method) {
    for (r in 1:embedding_dim) {
      W[, r, ] <- W_init[, r, ]
      lr <- step_size

      # parameters to judge the convergence.
      l <- 0
      record <- 0
      record_no <- 0
      ee <- 0

      # Initialize Adam parameters.
      if (update_method == "Adam") {
        momentum <- 0
        variance <- 0
        if (is.null(adam_parameter_init)) {
          beta1 <- 0.9
          beta2 <- 0.999
          epsilon_adam <- 1e-8
        } else {
          beta1 <- adam_parameter_init[1]
          beta2 <- adam_parameter_init[2]
          epsilon_adam <- adam_parameter_init[3]
        }
      }

      for (h in 1:iteration_step) {
        # loading the batch for this iteration step.
        batch_no <- ceiling(n_nodes / batch_size)
        if (h %% batch_no == 1 || batch_no == 1) {
          sample_index <- sample(n_nodes, n_nodes, replace = FALSE)
        }

        if (batch_no == 1) {
          mini_index_lower <- 1
          mini_index_upper <- n_nodes
        } else {
          mini_index_lower <- max(((h - 1) %% batch_no) * batch_size, 1)
          mini_index_upper <- min((h %% batch_no) * batch_size, n_nodes)
          if (mini_index_upper == 0) {
            mini_index_upper <- n_nodes
          }
          if (mini_index_upper - mini_index_lower == 1) {
            mini_index_lower <- mini_index_lower - round(batch_size / 2)
          }
        }
        mini_index <- sample_index[mini_index_lower:mini_index_upper]

        # stochastic gradient step.
        if (update_method == "Adam") {
          adam_parameter <- list(beta1, beta2, h, momentum, variance, epsilon_adam)
          adam_result <- gradient_computation(adjacency_tensor, W, B, BB, r, method = "Adam", mini_index = mini_index, adam_par = adam_parameter)
          W[, r, ] <- W[, r, ] - lr * adam_result[[1]]
          momentum <- adam_result[[2]]
          variance <- adam_result[[3]]
          rm(adam_result)
          gc()

          # calculate the objective and determine the convergence.
          if (h == 1) {
            l_old <- l
            l <- objective_computation(adjacency_tensor, W, B)
            l_record <- l
            if (!parallel) cat("iter=", h, "; ", "l=", l, "\n")
          } else if (h %% batch_no == 0) {
            # rm(derivate)
            gc()

            if (h %% (epoch_length * batch_no) == 0) {
              l_old <- l
              l <- objective_computation(adjacency_tensor, W, B)
              if (!parallel) cat("iter=", h, "; ", "l-l_old=", l - l_old, "\n")
              if (abs((l - l_old) / l_old) < epoch_length * epsilon) {
                if (!parallel) cat("final=", h, "; ", "l=", l, "\n")
                break
              }
            }

            if (h %% (10 * batch_no) == 0) {
              lr <- lr / 10
            }
          }
        } else {
          adam_result <- gradient_computation(adjacency_tensor, W, B, BB, r, method = "GD", mini_index = mini_index)
          W[, r, ] <- W[, r, ] - step_size * adam_result
          rm(adam_result)
          gc()

          # calculate the objective and determine the convergence.
          if (h == 1) {
            l_old <- l
            l <- objective_computation(adjacency_tensor, W, B)
            l_record <- l
            if (!parallel) cat("iter=", h, "; ", "l=", l, "\n")
          } else if (h %% (epoch_length) == 0) {
            l_old <- l
            l <- objective_computation(adjacency_tensor, W, B)
            if (!parallel) cat("iter=", h, "; ", "l-l_old=", l - l_old, "\n")
            if (abs((l - l_old) / l_old) < epoch_length * epsilon) {
              if (!parallel) cat("final=", h, "; ", "l=", l, "\n")
              break
            }
          }
        }
      }
    }
  } else if (!scalable_method) {
    batch_no <- 1
    for (r in 1:embedding_dim) {
      W[, r, ] <- W_init[, r, ]
      lr <- step_size

      # parameters to judge the convergence.
      l <- 0
      record <- 0
      record_no <- 0
      ee <- 0


      # Initialize Adam parameters.
      if (update_method == "Adam") {
        momentum <- 0
        variance <- 0
        if (is.null(adam_parameter_init)) {
          beta1 <- 0.9
          beta2 <- 0.999
          epsilon_adam <- 1e-8
        } else {
          beta1 <- adam_parameter_init[1]
          beta2 <- adam_parameter_init[2]
          epsilon_adam <- adam_parameter_init[3]
        }
      }

      for (h in 1:iteration_step) {
        # stochastic gradient step.
        if (update_method == "Adam") {
          adam_parameter <- list(beta1, beta2, h, momentum, variance, epsilon_adam)
          adam_result <- gradient_computation(adjacency_tensor, W, B, BB, r, method = "Adam", adam_par = adam_parameter)
          W[, r, ] <- W[, r, ] - lr * adam_result[[1]]
          momentum <- adam_result[[2]]
          variance <- adam_result[[3]]
          rm(adam_result)
          gc()

          # calculate the objective and determine the convergence.
          if (h == 1) {
            l_old <- l
            l <- objective_computation(adjacency_tensor, W, B)
            l_record <- l
            if (!parallel) cat("iter=", h, "; ", "l=", l, "\n")
          } else if (h %% batch_no == 0) {
            # rm(derivate)
            gc()

            if (h %% (epoch_length * batch_no) == 0) {
              l_old <- l
              l <- objective_computation(adjacency_tensor, W, B)
              if (!parallel) cat("iter=", h, "; ", "l-l_old=", l - l_old, "\n")
              if (abs((l - l_old) / l_old) < epoch_length * epsilon) {
                if (!parallel) cat("final=", h, "; ", "l=", l, "\n")
                break
              }
            }

            if (h %% (10 * batch_no) == 0) {
              lr <- lr / 10
            }
          }
        } else {
          adam_result <- gradient_computation(adjacency_tensor, W, B, BB, r, method = "GD")
          W[, r, ] <- W[, r, ] - step_size * adam_result
          rm(adam_result)
          gc()

          # calculate the objective and determine the convergence.
          if (h == 1) {
            l_old <- l
            l <- objective_computation(adjacency_tensor, W, B)
            l_record <- l
            if (!parallel) cat("iter=", h, "; ", "l=", l, "\n")
          } else {
            l_old <- l
            l <- objective_computation(adjacency_tensor, W, B)
            if (!parallel) cat("iter=", h, "; ", "l-l_old=", l - l_old, "\n")
            if (abs((l - l_old) / l_old) < epsilon) {
              if (!parallel) cat("final=", h, "; ", "l=", l, "\n")
              break
            }
          }
        }
      }
    }
  }

  W_T <- aperm(W, c(3, 1, 2))
  Z <- fda::fd(coef = W_T, basisobj = B_splines)
  list(Z, W, B, B_splines)
}

#' @title Perform (parallelized) hyperparameter tuning for network embeddings.
#'
#' @param spline_no_range Numeric vector, the range for tuning the number of B-splines.
#' @param embedding_dim_range Numeric vector, the range for tuning the dimension of latent functions.
#' @param tuning_method String, if tuning_method = "heuristic" then alternating optimizing the two hyperparameters, else if tuning_method = "grid" using the grid search.
#' @param adjacency_tensor an n*n*T array where each slice is a snapshot of the dynamic network.
#' @param spline_order  (optional) scalar. the polynomial order of B-spline basis, default 4
#' @param batch_size  (optional) scalar. the size of batches in mini-batch algorithm, default max(50,floor(n/20))
#' @param kernel_scale  (optional) scalar. the scaled parameter in the kernel used for smoothing, default two times the quartile of the time intervals.
#' @param timestamp_vec  (optional) vector. the observed time points of observation, default 1:T.
#' @param scalable_method  (optional) logical. whether to use the random-projection based SVD and mini-batch methods, default TRUE.
#' @param scalable_dim (optional) scalar. the dimension of random projection if scalable_method = TRUE, default 20.
#' @param scalable_power (optional) scalar. the time of power iterations in power method when using the random-projection based SVD, default 6.
#' @param iteration_step  (optional) scalar. the number of maximal iteration steps for the gradient algorithm, default 2500.
#' @param epsilon (optional) scalar. the stopping criterion for the gradient algorithm, default 1e-6.
#' @param parallel  (optional) logical. whether to parallel when tuning hyperparameters, defult TRUE.
#'
#' @return a matrix containing three columns, corresponding to spline numbers, embedding dimensions and the resulting generalized cross-validation criterion.
#' @export
#'
#' @examples
#' \donttest{
#' ER_generation <- function(n, p) {
#'   adjacency_matrix <- matrix(0, nrow = n, ncol = n)
#'   for (i in 2:n) {
#'     for (j in 1:(i - 1)) {
#'       adjacency_matrix[i, j] <- sample(c(0, 1), size = 1, prob = c(1 - p, p))
#'     }
#'     adjacency_matrix[i, i] <- 0.5
#'   }
#'   adjacency_matrix + t(adjacency_matrix)
#' }
#' set.seed(10)
#' n_nodes <- 10
#' T <- 10
#' p <- 0.2
#' dynamic_adjacency <- c()
#' dynamic_network_adjacency <- array(0, dim = c(n_nodes, n_nodes, T))
#' for (t in 1:T) {
#'   adj <- ER_generation(n_nodes, p)
#'   dynamic_adjacency <- c(dynamic_adjacency, list(list(adj)))
#'   dynamic_network_adjacency[, , t] <- adj
#' }
#' tuning_sim1 <- fASE_tuning(5:8, 4:5, tuning_method = "heuristic", adjacency_tensor = dynamic_network_adjacency, spline_order = 4, batch_size = NULL, kernel_scale = FALSE, timestamp_vec = NULL, scalable_method = TRUE, scalable_dim = 30, scalable_power = 6, iteration_step = 2500, epsilon = 1e-3)
#' tuning_parameter1 <- get_tuning(tuning_sim1)
#' }
fASE_tuning <- function(spline_no_range, embedding_dim_range, tuning_method = "heuristic", adjacency_tensor, spline_order = 4, batch_size = NULL, kernel_scale = FALSE, timestamp_vec = NULL, scalable_method = TRUE,
                        scalable_dim = 20, scalable_power = 6, iteration_step = 2000, epsilon = 1e-5, parallel = FALSE) {
  fASE <- function(adjacency_tensor, embedding_dim, spline_no, spline_order = 4, batch_size = NULL, kernel_scale = FALSE, timestamp_vec = NULL, scalable_method = TRUE,
                   scalable_dim = 20, scalable_power = 6, update_method = "Adam", adam_parameter = NULL, iteration_step = 2500, epsilon = 1e-6, step_size = 0.1, epoch_length = 1, parallel = FALSE) {
    set.seed(1)

    n_nodes <- dim(adjacency_tensor)[1]
    m_layers <- dim(adjacency_tensor)[3]
    adam_parameter_init <- adam_parameter

    if (is.null(timestamp_vec)) {
      timestamp_vec <- 1:m_layers
    }

    if (spline_no - spline_order + 1 > m_layers) {
      cat("the number of knots should be less than the number of layers!")
      break
    }

    if (scalable_method) {
      if (is.null(batch_size)) {
        if (n_nodes < 50) {
          batch_size <- n_nodes
        } else {
          batch_size <- max(50, floor(n_nodes / 20))
        }
      }
    }

    if (is.null(timestamp_vec)) {
      knots <- seq(1, m_layers, by = (m_layers - 1) / (spline_no - spline_order + 1))
    } else {
      knots <- seq(min(timestamp_vec), max(timestamp_vec), by = (max(timestamp_vec) - min(timestamp_vec)) / (spline_no - spline_order + 1))
    }

    # functions for initialization, calculation of objective functions and gradients.
    fx.hat <- function(x, z, h) {
      dnorm((z - x) / h) / h
    }
    smoothing <- function(y, a) {
      l <- length(a)
      new_mat <- matrix(0, nrow = dim(y)[1], ncol = dim(y)[2])
      for (i in 1:l) {
        new_mat <- new_mat + y[, , i] * a[i]
      }
      new_mat
    }
    KSMOOTH <- function(h, y, x, z) {
      m <- length(x)
      s.hat <- array(0, dim = c(dim(y)[1], dim(y)[2], m))
      for (i in 1:m) {
        a <- fx.hat(x[i], z, h)
        s.hat[, , i] <- smoothing(y, a)
      }
      s.hat
    }
    half_decompose <- function(X, d, s, q, scalable_method) {
      if (scalable_method) {
        set.seed(1)

        n <- dim(X)[1]

        Omega <- matrix(rnorm(n * s), nrow = n, ncol = s)
        A <- Omega

        for (i in 1:q) {
          if (i %% 2 == 1) {
            A <- X %*% A
          } else {
            A <- t(X) %*% A
          }
        }

        qr_result <- qr(A)
        rm(A)
        gc()
        Q <- qr.Q(qr_result)

        X_approx <- t(Q) %*% X %*% Q
        X_approx[which(X_approx != t(X_approx))] <- 0
        eigen_result <- eigen(X_approx)
        rm(X, X_approx)
        gc()

        eigenval <- eigen_result$values[1:d]
        eigenval[eigenval < 0] <- 0

        result <- t(t(Q %*% eigen_result$vectors[, 1:d]) * sqrt(eigenval))
      } else {
        eigen_result <- eigen(X)
        rm(X)
        gc()

        eigenval <- eigen_result$values[1:d]
        eigenval[eigenval < 0] <- 0

        result <- t(t(eigen_result$vectors[, 1:d]) * sqrt(eigenval))
      }

      result
    }
    kernel_embedding_initial <- function(adjacency_tensor, n_nodes, m_layers, timestamp_vec, embedding_dim, spline_no, B_splines, scalable_method, random_dim, mul_iter, h, kernel = "KS", spline_order = 4) {
      grid_no <- spline_no - spline_order + 2

      grids <- seq(min(timestamp_vec), max(timestamp_vec), by = (m_layers - 1) / (grid_no - 1))
      B <- t(fda::eval.basis(grids, B_splines))

      if (kernel == "KS") {
        if (!h) {
          h <- 2 * quantile(diff(timestamp_vec), 0.25)
        }
        K <- KSMOOTH(h, adjacency_tensor, grids, timestamp_vec)
      }
      if (is.null(random_dim)) {
        random_dim <- min(max(100, 10 * embedding_dim), n_nodes)
      }
      if (is.null(mul_iter)) {
        mul_iter <- 5
      }

      Z_init <- array(0, dim = c(n_nodes, embedding_dim, grid_no))

      for (l in 1:grid_no) {
        Z_init2 <- half_decompose(K[, , l], embedding_dim, random_dim, mul_iter, scalable_method)

        if (l > 1) {
          svd_result <- svd(t(Z_init2) %*% Z_init[, , l - 1])
          rotation_matrix <- svd_result$v %*% t(svd_result$u)
          Z_init[, , l] <- Z_init2 %*% rotation_matrix
        } else {
          Z_init[, , l] <- Z_init2
        }
      }

      W_init <- apply(Z_init, c(1, 2), crossprod, solve(t(B) %*% B, t(B)))
      W_init <- aperm(W_init, c(2, 3, 1))
      W_init
    }
    gradient_computation <- function(adjacency_tensor, W, B, BB, r, method = "GD", mini_index = NULL, adam_par = NULL) {
      m <- dim(adjacency_tensor)[3]
      n <- dim(adjacency_tensor)[1]
      d <- dim(W)[2]
      if (is.null(mini_index)) {
        mini_index <- 1:n
      }
      b <- length(mini_index)

      sum_matrices <- function(x, y) {
        x + y
      }


      if (method == "Adam") {
        # parameters of Adam are:
        # beta1, beta2, iter, m_(t-1), v_(t-1), epsilon
        beta1 <- adam_par[[1]]
        beta2 <- adam_par[[2]]
        iter <- adam_par[[3]]
        momentum <- adam_par[[4]]
        variance <- adam_par[[5]]
        epsilon <- adam_par[[6]]

        derivate <- 0
        WW <- array(apply(W, 2, crossprod, W[, r, ]), dim = c(dim(W)[3], dim(W)[3], dim(W)[2]))

        if (b <= 200) {
          for (i in 1:m) {
            temp <- adjacency_tensor[, , i]
            derivate <- derivate - n / b * temp[, mini_index] %*% ((W[, r, ] %*% BB[, , i])[mini_index, ])
          }
          for (i in 1:m) {
            for (s in 1:d) {
              derivate <- derivate + tcrossprod(W[, s, ], t(BB[, , i] %*% WW[, , s] %*% BB[, , i]))
            }
          }
        } else {
          for (i in 1:m) {
            temp <- adjacency_tensor[, , i]
            no_of_subbatches <- ceiling(b / 200)
            for (s in 1:no_of_subbatches) {
              mini_subindex <- (200 * (s - 1) + 1):min(200 * s, b)
              derivate <- derivate - n / b * temp[, mini_index[mini_subindex]] %*% ((W[, r, ] %*% BB[, , i])[mini_index[mini_subindex], ])
            }
          }
          for (i in 1:m) {
            for (s in 1:d) {
              derivate <- derivate + tcrossprod(W[, s, ], t(BB[, , i] %*% WW[, , s] %*% BB[, , i]))
            }
          }
        }


        # derivate = foreach(i=1:m, .combine = "sum_matrices")%dopar%{
        #   temp = adjacency_tensor[,,i]
        #   derivate = - n/b * temp[,mini_index]%*%((W[,r,]%*%B[,i]%*%t(B[,i]))[mini_index,])
        #   for(s in 1:d){
        #    derivate = derivate + W[,s,]%*%B[,i]%*%t(B[,i])%*%(t(W[,s,])%*%W[,r,])%*%B[,i]%*%t(B[,i])
        #   }
        #  derivate
        # }

        momentum <- beta1 * momentum + (1 - beta1) * derivate
        variance <- beta2 * variance + (1 - beta2) * derivate^2
        alpha_now <- sqrt(1 - beta2^iter) / (1 - beta1^iter)
        derivate <- alpha_now * momentum / (sqrt(variance) + epsilon)

        result <- list(derivate, momentum, variance)
      } else if (method == "GD") {
        derivate <- 0
        batch_size <- length(mini_index)

        WW <- array(apply(W, 2, crossprod, W[, r, ]), dim = c(dim(W)[3], dim(W)[3], dim(W)[2]))

        for (i in 1:m) {
          temp <- adjacency_tensor[, , i]
          no_of_batches <- ceiling(dim(adjacency_tensor)[1] / batch_size)
          for (s in 1:no_of_batches) {
            mini_index <- (batch_size * (s - 1) + 1):min(batch_size * s, dim(adjacency_tensor)[1])
            derivate <- derivate - temp[, mini_index] %*% ((W[, r, ] %*% B[, i] %*% t(B[, i]))[mini_index, ])
          }
          for (s in 1:d) {
            derivate <- derivate + W[, s, ] %*% B[, i] %*% t(B[, i]) %*% WW[, , s] %*% B[, i] %*% t(B[, i])
          }
        }

        result <- derivate
      }

      result
    }
    objective_computation <- function(adjacency_tensor, W, B) {
      objective_result <- 0

      m <- dim(adjacency_tensor)[3]
      n <- dim(adjacency_tensor)[1]
      d <- dim(W)[2]

      for (i in 1:m) {
        temp <- adjacency_tensor[, , i]
        for (s in 1:d) {
          tempp <- W[, s, ] %*% B[, i]
          temp <- temp - tempp %*% t(tempp)
          rm(tempp)
          gc()
        }

        objective_result <- objective_result + norm(temp, type = "F")^2
        rm(temp)
        gc()
      }

      objective_result
    }


    B_splines <- fda::create.bspline.basis(rangeval = c(knots[1], knots[length(knots)]), nbasis = spline_no, norder = spline_order, breaks = knots)
    B <- t(fda::eval.basis(timestamp_vec, B_splines))
    BB <- array(apply(B, 2, tcrossprod), dim = c(dim(B)[1], dim(B)[1], dim(B)[2]))

    W_init <- kernel_embedding_initial(adjacency_tensor, n_nodes, m_layers, timestamp_vec, embedding_dim, spline_no, B_splines, scalable_method, scalable_dim, scalable_power, kernel_scale, spline_order = spline_order)
    W <- W_init

    if (scalable_method) {
      for (r in 1:embedding_dim) {
        W[, r, ] <- W_init[, r, ]
        lr <- step_size

        # parameters to judge the convergence.
        l <- 0
        record <- 0
        record_no <- 0
        ee <- 0

        # Initialize Adam parameters.
        if (update_method == "Adam") {
          momentum <- 0
          variance <- 0
          if (is.null(adam_parameter_init)) {
            beta1 <- 0.9
            beta2 <- 0.999
            epsilon_adam <- 1e-8
          } else {
            beta1 <- adam_parameter_init[1]
            beta2 <- adam_parameter_init[2]
            epsilon_adam <- adam_parameter_init[3]
          }
        }

        for (h in 1:iteration_step) {
          # loading the batch for this iteration step.
          batch_no <- ceiling(n_nodes / batch_size)
          if (h %% batch_no == 1 || batch_no == 1) {
            sample_index <- sample(n_nodes, n_nodes, replace = FALSE)
          }

          if (batch_no == 1) {
            mini_index_lower <- 1
            mini_index_upper <- n_nodes
          } else {
            mini_index_lower <- max(((h - 1) %% batch_no) * batch_size, 1)
            mini_index_upper <- min((h %% batch_no) * batch_size, n_nodes)
            if (mini_index_upper == 0) {
              mini_index_upper <- n_nodes
            }
            if (mini_index_upper - mini_index_lower == 1) {
              mini_index_lower <- mini_index_lower - round(batch_size / 2)
            }
          }
          mini_index <- sample_index[mini_index_lower:mini_index_upper]

          # stochastic gradient step.
          if (update_method == "Adam") {
            adam_parameter <- list(beta1, beta2, h, momentum, variance, epsilon_adam)
            adam_result <- gradient_computation(adjacency_tensor, W, B, BB, r, method = "Adam", mini_index = mini_index, adam_par = adam_parameter)
            W[, r, ] <- W[, r, ] - lr * adam_result[[1]]
            momentum <- adam_result[[2]]
            variance <- adam_result[[3]]
            rm(adam_result)
            gc()

            # calculate the objective and determine the convergence.
            if (h == 1) {
              l_old <- l
              l <- objective_computation(adjacency_tensor, W, B)
              l_record <- l
              if (!parallel) cat("iter=", h, "; ", "l=", l, "\n")
            } else if (h %% batch_no == 0) {
              # rm(derivate)
              gc()

              if (h %% (epoch_length * batch_no) == 0) {
                l_old <- l
                l <- objective_computation(adjacency_tensor, W, B)
                if (!parallel) cat("iter=", h, "; ", "l-l_old=", l - l_old, "\n")
                if (abs((l - l_old) / l_old) < epoch_length * epsilon) {
                  if (!parallel) cat("final=", h, "; ", "l=", l, "\n")
                  break
                }
              }

              if (h %% (10 * batch_no) == 0) {
                lr <- lr / 10
              }
            }
          } else {
            adam_result <- gradient_computation(adjacency_tensor, W, B, BB, r, method = "GD", mini_index = mini_index)
            W[, r, ] <- W[, r, ] - step_size * adam_result
            rm(adam_result)
            gc()

            # calculate the objective and determine the convergence.
            if (h == 1) {
              l_old <- l
              l <- objective_computation(adjacency_tensor, W, B)
              l_record <- l
              if (!parallel) cat("iter=", h, "; ", "l=", l, "\n")
            } else if (h %% (epoch_length) == 0) {
              l_old <- l
              l <- objective_computation(adjacency_tensor, W, B)
              if (!parallel) cat("iter=", h, "; ", "l-l_old=", l - l_old, "\n")
              if (abs((l - l_old) / l_old) < epoch_length * epsilon) {
                if (!parallel) cat("final=", h, "; ", "l=", l, "\n")
                break
              }
            }
          }
        }
      }
    } else if (!scalable_method) {
      batch_no <- 1
      for (r in 1:embedding_dim) {
        W[, r, ] <- W_init[, r, ]
        lr <- step_size

        # parameters to judge the convergence.
        l <- 0
        record <- 0
        record_no <- 0
        ee <- 0


        # Initialize Adam parameters.
        if (update_method == "Adam") {
          momentum <- 0
          variance <- 0
          if (is.null(adam_parameter_init)) {
            beta1 <- 0.9
            beta2 <- 0.999
            epsilon_adam <- 1e-8
          } else {
            beta1 <- adam_parameter_init[1]
            beta2 <- adam_parameter_init[2]
            epsilon_adam <- adam_parameter_init[3]
          }
        }

        for (h in 1:iteration_step) {
          # stochastic gradient step.
          if (update_method == "Adam") {
            adam_parameter <- list(beta1, beta2, h, momentum, variance, epsilon_adam)
            adam_result <- gradient_computation(adjacency_tensor, W, B, BB, r, method = "Adam", adam_par = adam_parameter)
            W[, r, ] <- W[, r, ] - lr * adam_result[[1]]
            momentum <- adam_result[[2]]
            variance <- adam_result[[3]]
            rm(adam_result)
            gc()

            # calculate the objective and determine the convergence.
            if (h == 1) {
              l_old <- l
              l <- objective_computation(adjacency_tensor, W, B)
              l_record <- l
              if (!parallel) cat("iter=", h, "; ", "l=", l, "\n")
            } else if (h %% batch_no == 0) {
              # rm(derivate)
              gc()

              if (h %% (epoch_length * batch_no) == 0) {
                l_old <- l
                l <- objective_computation(adjacency_tensor, W, B)
                if (!parallel) cat("iter=", h, "; ", "l-l_old=", l - l_old, "\n")
                if (abs((l - l_old) / l_old) < epoch_length * epsilon) {
                  if (!parallel) cat("final=", h, "; ", "l=", l, "\n")
                  break
                }
              }

              if (h %% (10 * batch_no) == 0) {
                lr <- lr / 10
              }
            }
          } else {
            adam_result <- gradient_computation(adjacency_tensor, W, B, BB, r, method = "GD")
            W[, r, ] <- W[, r, ] - step_size * adam_result
            rm(adam_result)
            gc()

            # calculate the objective and determine the convergence.
            if (h == 1) {
              l_old <- l
              l <- objective_computation(adjacency_tensor, W, B)
              l_record <- l
              if (!parallel) cat("iter=", h, "; ", "l=", l, "\n")
            } else {
              l_old <- l
              l <- objective_computation(adjacency_tensor, W, B)
              if (!parallel) cat("iter=", h, "; ", "l-l_old=", l - l_old, "\n")
              if (abs((l - l_old) / l_old) < epsilon) {
                if (!parallel) cat("final=", h, "; ", "l=", l, "\n")
                break
              }
            }
          }
        }
      }
    }

    W_T <- aperm(W, c(3, 1, 2))
    Z <- fda::fd(coef = W_T, basisobj = B_splines)
    list(Z, W, B, B_splines)
  }
  objective_computation <- function(adjacency_tensor, W, B) {
    objective_result <- 0

    m <- dim(adjacency_tensor)[3]
    n <- dim(adjacency_tensor)[1]
    d <- dim(W)[2]

    for (i in 1:m) {
      temp <- adjacency_tensor[, , i]
      for (s in 1:d) {
        tempp <- W[, s, ] %*% B[, i]
        temp <- temp - tempp %*% t(tempp)
        rm(tempp)
        gc()
      }

      objective_result <- objective_result + norm(temp, type = "F")^2
      rm(temp)
      gc()
    }

    objective_result
  }

  # tuning important parameters spline_no, embedding_dim.
  n_nodes <- dim(adjacency_tensor)[1]
  m_layers <- dim(adjacency_tensor)[3]
  NGCV_list <- c()

  if (tuning_method == "grid") {
    if (!parallel) {
      for (spline_no in spline_no_range) {
        for (embedding_dim in embedding_dim_range) {
          result <- fASE(adjacency_tensor, embedding_dim, spline_no, spline_order = spline_order, batch_size = batch_size, timestamp_vec = timestamp_vec, scalable_dim = scalable_dim, kernel_scale = kernel_scale, iteration_step = iteration_step, epsilon = epsilon)
          W <- result[[2]]
          B <- result[[3]]

          log_l <- log(objective_computation(adjacency_tensor, W, B) / (n_nodes^2 * m_layers))
          NGCV <- log_l - 2 * log(1 - (2 * spline_no * embedding_dim) / (n_nodes * m_layers))

          NGCV_list <- rbind(NGCV_list, t(c(spline_no, embedding_dim, NGCV)))

          cat("(", spline_no, ",", embedding_dim, ") finished.\n")

          # save(file = "result.RData", NGCV_list)
        }
      }
    } else {
      numCores <- parallel::detectCores()
      cl <- makeCluster(numCores)
      registerDoParallel(cl)

      NGCV_list <- foreach(spline_no = spline_no_range, .combine = "rbind", .verbose = TRUE) %:%
        foreach(embedding_dim = embedding_dim_range, .combine = "rbind", .verbose = TRUE) %dopar% {
          result <- fASE(adjacency_tensor, embedding_dim, spline_no, spline_order = spline_order, batch_size = batch_size, timestamp_vec = timestamp_vec, scalable_dim = scalable_dim, kernel_scale = kernel_scale, iteration_step = iteration_step, epsilon = epsilon, parallel = TRUE)
          W <- result[[2]]
          B <- result[[3]]

          log_l <- log(objective_computation(adjacency_tensor, W, B) / (n_nodes^2 * m_layers))
          NGCV <- log_l - 2 * log(1 - (2 * spline_no * embedding_dim) / (n_nodes * m_layers))

          cat("(", spline_no, ",", embedding_dim, ") finished.\n")

          t(c(spline_no, embedding_dim, NGCV))

          # save(file = "result.RData", NGCV_list)
        }
      stopCluster(cl)
    }
  } else {
    embedding_dim <- embedding_dim_range[1]
    old_min_NGCV <- Inf
    min_NGCV <- Inf

    if (!parallel) {
      while (old_min_NGCV > min_NGCV || is.infinite(min_NGCV)) {
        old_min_NGCV <- min_NGCV

        for (spline_no in spline_no_range) {
          result <- fASE(adjacency_tensor, embedding_dim, spline_no, spline_order = spline_order, batch_size = batch_size, timestamp_vec = timestamp_vec, scalable_dim = scalable_dim, kernel_scale = kernel_scale, iteration_step = iteration_step, epsilon = epsilon)
          W <- result[[2]]
          B <- result[[3]]

          log_l <- log(objective_computation(adjacency_tensor, W, B) / (n_nodes^2 * m_layers))
          NGCV <- log_l - 2 * log(1 - (2 * spline_no * embedding_dim) / (n_nodes * m_layers))

          if (NGCV < min_NGCV) {
            min_NGCV <- NGCV
          }

          NGCV_list <- rbind(NGCV_list, t(c(spline_no, embedding_dim, NGCV)))

          cat("(", spline_no, ",", embedding_dim, ") finished.\n")
        }

        spline_no <- NGCV_list[which(NGCV_list[, 3] == min_NGCV)[1], 1]

        if (min_NGCV == old_min_NGCV) {
          break
        }

        for (embedding_dim in embedding_dim_range) {
          result <- fASE(adjacency_tensor, embedding_dim, spline_no, spline_order = spline_order, batch_size = batch_size, timestamp_vec = timestamp_vec, scalable_dim = scalable_dim, kernel_scale = kernel_scale, iteration_step = iteration_step, epsilon = epsilon)
          W <- result[[2]]
          B <- result[[3]]

          log_l <- log(objective_computation(adjacency_tensor, W, B) / (n_nodes^2 * m_layers))
          NGCV <- log_l - 2 * log(1 - (2 * spline_no * embedding_dim) / (n_nodes * m_layers))

          if (NGCV < min_NGCV) {
            min_NGCV <- NGCV
          }

          NGCV_list <- rbind(NGCV_list, t(c(spline_no, embedding_dim, NGCV)))

          cat("(", spline_no, ",", embedding_dim, ") finished.\n")
        }

        embedding_dim <- NGCV_list[which(NGCV_list[, 3] == min_NGCV)[1], 2]
      }
    } else {
      numCores <- parallel::detectCores()
      cl <- makeCluster(numCores)
      registerDoParallel(cl)

      while (old_min_NGCV > min_NGCV || is.infinite(min_NGCV)) {
        old_min_NGCV <- min_NGCV

        NGCV_list2 <- foreach(spline_no = spline_no_range, .combine = "rbind", .verbose = TRUE) %dopar% {
          result <- fASE(adjacency_tensor, embedding_dim, spline_no, spline_order = spline_order, batch_size = batch_size, timestamp_vec = timestamp_vec, scalable_dim = scalable_dim, kernel_scale = kernel_scale, iteration_step = iteration_step, epsilon = epsilon, parallel = TRUE)
          W <- result[[2]]
          B <- result[[3]]

          log_l <- log(objective_computation(adjacency_tensor, W, B) / (n_nodes^2 * m_layers))
          NGCV <- log_l - 2 * log(1 - (2 * spline_no * embedding_dim) / (n_nodes * m_layers))

          cat("(", spline_no, ",", embedding_dim, ") finished.\n")

          t(c(spline_no, embedding_dim, NGCV))
        }

        min_NGCV <- min(min(NGCV_list2[, 3]), min_NGCV)
        NGCV_list <- rbind(NGCV_list, NGCV_list2)

        spline_no <- NGCV_list[which(NGCV_list[, 3] == min_NGCV)[1], 1]

        if (min_NGCV == old_min_NGCV) {
          break
        }

        NGCV_list2 <- foreach(embedding_dim = embedding_dim_range, .combine = "rbind", .verbose = TRUE) %dopar% {
          result <- fASE(adjacency_tensor, embedding_dim, spline_no, spline_order = spline_order, batch_size = batch_size, timestamp_vec = timestamp_vec, scalable_dim = scalable_dim, kernel_scale = kernel_scale, iteration_step = iteration_step, epsilon = epsilon, parallel = TRUE)
          W <- result[[2]]
          B <- result[[3]]

          log_l <- log(objective_computation(adjacency_tensor, W, B) / (n_nodes^2 * m_layers))
          NGCV <- log_l - 2 * log(1 - (2 * spline_no * embedding_dim) / (n_nodes * m_layers))

          cat("(", spline_no, ",", embedding_dim, ") finished.\n")

          t(c(spline_no, embedding_dim, NGCV))
        }

        min_NGCV <- min(min(NGCV_list2[, 3]), min_NGCV)
        NGCV_list <- rbind(NGCV_list, NGCV_list2)

        embedding_dim <- NGCV_list[which(NGCV_list[, 3] == min_NGCV)[1], 2]
      }
      stopCluster(cl)
    }
  }

  NGCV_list
}

#' @title Get the optimal hyperparameters from an NGCV list.
#'
#' @param NGCV_list the matrix obtained by the fASE_tuning function.
#'
#' @return a list of optimal hyperparameters
#' @export
#'
#' @examples
#' \donttest{
#' ER_generation <- function(n, p) {
#'   adjacency_matrix <- matrix(0, nrow = n, ncol = n)
#'   for (i in 2:n) {
#'     for (j in 1:(i - 1)) {
#'       adjacency_matrix[i, j] <- sample(c(0, 1), size = 1, prob = c(1 - p, p))
#'     }
#'     adjacency_matrix[i, i] <- 0.5
#'   }
#'   adjacency_matrix + t(adjacency_matrix)
#' }
#' set.seed(10)
#' n_nodes <- 10
#' T <- 10
#' p <- 0.2
#' dynamic_adjacency <- c()
#' dynamic_network_adjacency <- array(0, dim = c(n_nodes, n_nodes, T))
#' for (t in 1:T) {
#'   adj <- ER_generation(n_nodes, p)
#'   dynamic_adjacency <- c(dynamic_adjacency, list(list(adj)))
#'   dynamic_network_adjacency[, , t] <- adj
#' }
#' tuning_sim1 <- fASE_tuning(5:8, 4:5, tuning_method = "heuristic", adjacency_tensor = dynamic_network_adjacency, spline_order = 4, batch_size = NULL, kernel_scale = FALSE, timestamp_vec = NULL, scalable_method = TRUE, scalable_dim = 30, scalable_power = 6, iteration_step = 2500, epsilon = 1e-3)
#' tuning_parameter1 <- get_tuning(tuning_sim1)
#' }
get_tuning <- function(NGCV_list) {
  index <- which(NGCV_list[, 3] == min(NGCV_list[, 3]))[1]
  list(NGCV_list[index, 1], NGCV_list[index, 2], NGCV_list[index, 3])
}
