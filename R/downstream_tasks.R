# library(igraph)
# #library(fda)
# library(fda.usc)
# library(MASS)
# #library(ggplot2)
# #library(ggalt)
# library(grid)
# #library(reshape2)
# library(patchwork)
# library(lubridate)
# #library(factoextra)
# library(fpc)
# #library(dplyr)
# library(tidyr)
# library(readr)
# # library(MFPCA)
# library(roahd)
# library(ecp)
# library(cluster)
# #library(kselection)
# #library(lowmemtkmeans)
# #library(doSNOW)
# library(robustbase)
# library(Rtsne)
library(dplyr)
library(ggplot2)

#' @title Conduct functional PCA for latent functions without radius correction
#' @param embedding_fASE a T*n*d array. values of the resulting latent functions obtained by fASE function at timestamp_vec.
#' @param timestamp_vec T vector. the observed time points of observation.
#' @param PCA_dim scalar. the functional PCA dimension for latent functions.
#' @param active_points (optional) list. If the option is not NULL, it must be a list containing active time points of each node.
#'
#' @return a list of three elements: fPC scores, functional principal components, proportion of variations explained
#' @export
#'
#' @examples
#' data("SBM_example_fASE0.rda")
#' result_SBM0 <- SBM_example_fASE0[[1]]
#' embedding_SBM0 <- fda::eval.fd(1:50, result_SBM0[[1]])
#' PCA_score <- fPCA(embedding_SBM0, 1:50, 2)
fPCA <- function(embedding_fASE, timestamp_vec, PCA_dim, active_points = NULL) {
  m_layers <- dim(embedding_fASE)[1]
  embedding_dim <- dim(embedding_fASE)[3]
  if (!is.null(active_points)) {
    active_number <- sapply(1:dim(embedding_fASE)[2], function(x) {
      length(active_points[[x]])
    })
    active_number[active_number == 0] <- 1
    embedding_fASE <- embedding_fASE / sqrt(aperm(array(active_number, dim = c(dim(embedding_fASE)[2], dim(embedding_fASE)[1], dim(embedding_fASE)[3])), c(2, 1, 3)))
  }

  # the third dim of embedding_fASE should be consistent with the length of timestamp_vec.

  if (length(timestamp_vec) != m_layers) {
    cat("Parameters do not fit with each other!\n")
    break
  } else {
    ZZ <- vector(mode = "list", length = embedding_dim)
    uniExpansions <- vector(mode = "list", length = embedding_dim)
    for (i in 1:embedding_dim) {
      if (is.matrix(embedding_fASE[, , i])) {
        Zi <- funData::funData(argvals = timestamp_vec, X = t(embedding_fASE[, , i]))
      } else {
        Zi <- funData::funData(argvals = timestamp_vec, X = as.matrix(embedding_fASE[, , i]))
      }

      ZZ[[i]] <- Zi
      uniExpansions[[i]] <- list(type = "fda")
    }

    ZZ <- funData::multiFunData(ZZ)
    # pca_fASE = MFPCA::MFPCA(ZZ, PCA_dim, uniExpansions, weights = cumprod(rep(0.8,times=length(ZZ)))/sum(cumprod(rep(0.8,times=length(ZZ)))))
    pca_fASE <- MFPCA::MFPCA(ZZ, PCA_dim, uniExpansions)

    list(pca_fASE$scores, pca_fASE$functions, summary(pca_fASE)["Proportion of variance explained", ])
  }
}


#' @title Conduct functional PCA for latent functions with radius correction
#' @param embedding_fASE an T*n*d array. values of the resulting latent functions obtained by fASE function at timestamp_vec.
#' @param timestamp_vec T vector. the observed time points of observation.
#' @param PCA_dim scalar. the functional PCA dimension for latent functions.
#' @param active_points (optional) list. If the option is not NULL, it must be a list containing active time points of each node.
#'
#' @return a list of three elements: fPC scores, functional principal components, proportion of variations explained
#' @export
#'
#' @examples
#' data("DCBM_example_fASE0.rda")
#' result_DCBM0 <- DCBM_example_fASE0[[1]]
#' deg <- DCBM_example_fASE0[[2]]
#' embedding_DCBM0 <- fda::eval.fd(1:100, result_DCBM0[[1]])
#' PCA_score <- fPCA_angle(embedding_DCBM0, 1:100, 2)
fPCA_angle <- function(embedding_fASE, timestamp_vec, PCA_dim, active_points = NULL) {
  sintegral_use <- function(fx, x) {
    Bolstad2::sintegral(x, fx)$int
  }
  subset <- function(i, ls, x) {
    if (length(ls[[i]]) != 0) {
      h <- matrix(0, nrow = dim(x)[1], ncol = dim(x)[3])
      h[which((1:length(timestamp_vec)) %in% ls[[i]]), ] <- x[which((1:length(timestamp_vec)) %in% ls[[i]]), i, ]
    } else {
      h <- matrix(0, nrow = dim(x)[1], ncol = dim(x)[3])
    }
    h
  }

  embedding_fASE2 <- embedding_fASE
  if (!is.null(active_points)) {
    active_points <- sapply(1:length(active_points), function(i) {
      pmatch(active_points[[i]], timestamp_vec)
    })
    embedding_fASE2 <- sapply(1:dim(embedding_fASE)[2], subset, active_points, embedding_fASE2, simplify = "array")
    embedding_fASE2 <- aperm(embedding_fASE2, c(1, 3, 2))
  }
  embedding_fASE3 <- apply(embedding_fASE2, c(1, 2), crossprod)
  embedding_fASE3[embedding_fASE3 == 0] <- 1

  embedding_fASE <- embedding_fASE2 / sqrt(array(embedding_fASE3, dim = c(dim(embedding_fASE)[1], dim(embedding_fASE)[2], dim(embedding_fASE)[3])))
  if (!is.null(active_points)) {
    active_number <- sapply(1:dim(embedding_fASE)[2], function(x) {
      length(active_points[[x]])
    })
    active_number[active_number == 0] <- 1
    embedding_fASE <- embedding_fASE / sqrt(aperm(array(active_number, dim = c(dim(embedding_fASE)[2], dim(embedding_fASE)[1], dim(embedding_fASE)[3])), c(2, 1, 3)))
  }

  m_layers <- dim(embedding_fASE)[1]
  embedding_dim <- dim(embedding_fASE)[3]
  # the third dim of embedding_fASE should be consistent with the length of timestamp_vec.

  if (length(timestamp_vec) != m_layers) {
    cat("Parameters do not fit with each other!\n")
    break
  } else {
    ZZ <- vector(mode = "list", length = embedding_dim)
    uniExpansions <- vector(mode = "list", length = embedding_dim)
    for (i in 1:embedding_dim) {
      if (is.matrix(embedding_fASE[, , i])) {
        Zi <- funData::funData(argvals = timestamp_vec, X = t(embedding_fASE[, , i]))
      } else {
        Zi <- funData::funData(argvals = timestamp_vec, X = as.matrix(embedding_fASE[, , i]))
      }

      ZZ[[i]] <- Zi
      uniExpansions[[i]] <- list(type = "fda")
    }

    ZZ <- funData::multiFunData(ZZ)
    pca_fASE <- MFPCA::MFPCA(ZZ, PCA_dim, uniExpansions)

    list(pca_fASE$scores, pca_fASE$functions, summary(pca_fASE)["Proportion of variance explained", ])
  }
}

#' @title Calculate the active time points for each node
#' @description If the node has any links with others, we regard it as 'active'; otherwise, we regard it as 'inactive'.
#' @param node the node for the active time points calculation.
#' @param net an array where each slice is a snapshot of the dynamic network.
#' @param timestamp_vec (optional) the observed time points of observation.
#'
#' @return If timestamp_vec = NULL, output the indecies of snapshots where the node is active; otherwise, output the time points where the node is active
#' @export
#'
#' @examples
#' ER_varynodes <- function(n, p, t) {
#'   adjacency_matrix <- matrix(0, nrow = n, ncol = n)
#'
#'   active_nodes <- unique(c((floor((t - 1) / 10) * 10 + 1):((floor((t - 1) / 10) + 1) * 10), 80:100))
#'   for (i in active_nodes) {
#'     for (j in active_nodes) {
#'       if (i > j) {
#'         adjacency_matrix[i, j] <- sample(c(0, 1), size = 1, prob = c(1 - p, p))
#'       }
#'     }
#'     adjacency_matrix[i, i] <- 0.5
#'   }
#'   adjacency_matrix + t(adjacency_matrix)
#' }
#'
#' set.seed(10)
#' n_nodes <- 100
#' T <- 100
#' p <- 0.1
#' C <- c(rep(2, times = 79), rep(1, times = 21))
#' dynamic_adjacency <- c()
#' dynamic_network_adjacency <- array(0, dim = c(n_nodes, n_nodes, T))
#' for (t in 1:T) {
#'   adj <- ER_varynodes(n_nodes, p, t)
#'   dynamic_adjacency <- c(dynamic_adjacency, list(list(adj)))
#'   dynamic_network_adjacency[, , t] <- adj
#' }
#' active_points <- lapply(1:100, active_calculation, dynamic_network_adjacency)
active_calculation <- function(node, net, timestamp_vec = NULL) {
  if (is.null(timestamp_vec)) {
    which(apply(net[node, , ], 2, sum) != 0)
  } else {
    timestamp_vec[which(apply(net[node, , ], 2, sum) != 0)]
  }
}

#' @title Calculate cosine similarity
#' @description Calculate cosine similarity of latent functions for each pair of nodes in the entire time interval.
#' @param embedding_fASE a T*n*d array. values of the resulting latent functions obtained by fASE function at timestamp_vec.
#' @param timevec T vector. the observed time points of observation.
#'
#' @return an n*n matrix where each entry contains the cosine similarity of the corresponding node pair.
#' @export
#'
#' @examples
#' data("DCBM_example_fASE0.rda")
#' result_DCBM0 <- DCBM_example_fASE0[[1]]
#' deg <- DCBM_example_fASE0[[2]]
#' sim <- cosine_similarity(fda::eval.fd(seq(1, 100, length = 1000), result_DCBM0[[1]]), seq(1, 100, length = 1000))
cosine_similarity <- function(embedding_fASE, timevec) {
  sintegral_use <- function(fx, x) {
    Bolstad2::sintegral(x, fx)$int
  }
  cross_mul <- function(i, arr1, arr2) {
    res <- arr1[i, , ] %*% t(arr2[i, , ])
  }
  similarity_int <- mapply(cross_mul, i = 1:dim(embedding_fASE)[1], MoreArgs = list(arr1 = embedding_fASE, arr2 = embedding_fASE))
  similarity_mat <- apply(similarity_int, 1, sintegral_use, timevec)
  similarity_mat <- matrix(similarity_mat, nrow = dim(embedding_fASE)[2])
  similarity_mat <- similarity_mat / sqrt(diag(similarity_mat) %*% t(diag(similarity_mat)))
  similarity_mat
}

#' @title Calculate cosine similarity for all nodes in a batch
#' @description Calculate cosine similarity of latent functions for pre-assigned nodes and all other nodes, in several pre-assigned, different intervals.
#' @param embedding_fASE a T*n*d array. values of the resulting latent functions obtained by fASE function at timestamp_vec.
#' @param nodes_calc vector. indicies of nodes included in the calculation of the cosine similarity matrix.
#' @param timevec_list (optional if timevec is not NULL) list. contains several index vectors of time points in timevec, each of which are used as the time period for calculation of the cosine similarity matrix.
#' @param timevec (optional) T vector. the observed time points of observation.
#' @return a list of the same length as timevec_list, containing n*n  cosine similarity matrices for the specific time period.
#' @export
#'
#' @examples
#' data("DCBM_example_fASE0.RData")
#' result_DCBM0 <- DCBM_example_fASE0[[1]]
#' deg <- DCBM_example_fASE0[[2]]
#' sim <- cosine_similarity_partial(fda::eval.fd(c(seq(1, 50, length = 500), seq(51, 100, length = 500)), result_DCBM0[[1]]), 1:10, timevec_list = list(1:500, 501:1000))
cosine_similarity_partial <- function(embedding_fASE, nodes_calc, timevec_list = NULL, timevec = NULL) {
  if (is.null(timevec_list)) {
    if (!is.null(timevec)) {
      timevec_list <- as.list(timevec)
    } else {
      timevec_list <- as.list(1:dim(embedding_fASE)[1])
    }
  }
  sintegral_use <- function(fx, x) {
    if (length(x) > 1) {
      Bolstad2::sintegral(x, fx)$int
    } else {
      fx
    }
  }
  cross_mul <- function(i, arr1, arr2) {
    res <- arr1[i, , ] %*% t(arr2[i, , ])
  }
  similarity_mat <- list()
  for (s in 1:length(timevec_list)) {
    if (length(timevec_list[[s]]) != 1) {
      arr1 <- array(embedding_fASE[timevec_list[[s]], nodes_calc, ], dim = c(length(timevec_list[[s]]), length(nodes_calc), dim(embedding_fASE)[3]))
      arr2 <- embedding_fASE[timevec_list[[s]], , ]
    } else {
      arr1 <- aperm(array(embedding_fASE[timevec_list[[s]], nodes_calc, ], dim = c(length(nodes_calc), dim(embedding_fASE)[3], length(timevec_list[[s]]))), c(3, 1, 2))
      arr2 <- aperm(array(embedding_fASE[timevec_list[[s]], , ], dim = c(dim(embedding_fASE)[2], dim(embedding_fASE)[3], length(timevec_list[[s]]))), c(3, 1, 2))
    }
    similarity_int <- mapply(cross_mul, i = 1:dim(arr1)[1], MoreArgs = list(arr1 = arr1, arr2 = arr2))
    similarity_int <- array(similarity_int, dim = c(dim(arr1)[1], dim(arr1)[2], dim(arr2)[2]))
    norm_int <- apply(arr2, c(1, 2), function(x) {
      sum(x^2)
    })
    similarity_mat[[s]] <- apply(similarity_int, c(2, 3), sintegral_use, timevec_list[[s]])
    norm_vec <- apply(norm_int, 2, sintegral_use, timevec_list[[s]])
    similarity_mat[[s]] <- matrix(similarity_mat[[s]], nrow = length(nodes_calc))
    similarity_mat[[s]] <- similarity_mat[[s]] / sqrt(norm_vec[nodes_calc] %*% t(norm_vec))
  }

  similarity_mat
}


#' @title Calculate the functional depth using the cosine similarity matrix
#' @description Calculate the functional depth of each node using the cosine similarity matrix.
#' @param dist_mat an n*n matrix. the sine values of latent functions for each pair of nodes (1 - cosine_simlarity_matrix).
#'
#' @return vector. calculated depth for each node.
#' @export
#'
#' @examples
#' data("DCBM_example_fASE0.rda")
#' result_DCBM0 <- DCBM_example_fASE0[[1]]
#' deg <- DCBM_example_fASE0[[2]]
#' sim <- cosine_similarity(fda::eval.fd(1:100, result_DCBM0[[1]]), 1:100)
#' depth_stats <- metric_depth(1 - sim)
metric_depth <- function(dist_mat) {
  n_nodes <- dim(dist_mat)[1]
  depth_record <- rep(0, times = n_nodes)

  for (a in 1:n_nodes) {
    count_mat <- pmax(matrix(dist_mat[, a], nrow = dim(dist_mat)[1], ncol = dim(dist_mat)[2]), matrix(dist_mat[a, ], nrow = dim(dist_mat)[2], ncol = dim(dist_mat)[1], byrow = TRUE))
    count <- sum(dist_mat > count_mat) - sum(diag(dist_mat > count_mat))
    depth_record[a] <- count / (choose(n_nodes, 2) * 2)
  }
  depth_record
}

spherical_transform <- function(X) {
  m <- length(X)
  r <- sqrt(sum(X^2))
  phi <- rep(0, times = m)
  phi[1] <- r
  temp <- 1
  for (i in rev(1:(m - 1))) {
    phi[i + 1] <- atan(X[i + 1] / (X[i] * temp))
    temp <- temp * cos(phi[i + 1])
  }
  phi
}
radius_transform <- function(X) {
  Y <- apply(X, c(1, 2), function(x) {
    sqrt(sum(x^2))
  })
  Y[Y == 0] <- 1
  X / array(Y, dim = dim(X))
}

#' @title Perform K-means algorithm.
#' @description Perform K-means algorithm for a dataframe.
#' @param X an n*d data matrix.
#' @param k the number of clusters K.
#'
#' @return a list containing: cluster membership and withinss.
cluster_kmeans <- function(X, k) {
  set.seed(1)

  calculate_wcss <- function(X, cluster_centers, cluster_assignments) {
    wcss <- 0
    for (i in 1:nrow(X)) {
      # 获取数据点i所属的簇中心和簇分配
      cluster_center <- cluster_centers[cluster_assignments[i], ]
      # 计算欧氏距离的平方
      distance_squared <- sum((X[i, ] - cluster_center)^2)
      wcss <- wcss + distance_squared
    }
    return(wcss)
  }
  Kmeans_result <- ClusterR::KMeans_rcpp(X, k)
  clusters <- Kmeans_result$clusters
  withinss <- calculate_wcss(X, Kmeans_result$centroids, clusters)

  list(cluster = clusters, withinss = withinss)
}

#' @title Perform trimmed K-means algorithm
#' @description Perform trimmed K-means algorithm for a dataframe.
#' @param X an n*d data matrix.
#' @param k the number of clusters K.
#'
#' @return a list containing: cluster membership and withinss.
cluster_tkmeans <- function(X, k) {
  set.seed(1)
  calculate_wcss <- function(X, cluster_centers, cluster_assignments) {
    wcss <- 0
    for (i in 1:nrow(X)) {
      # 获取数据点i所属的簇中心和簇分配
      cluster_center <- cluster_centers[cluster_assignments[i], ]
      # 计算欧氏距离的平方
      distance_squared <- sum((X[i, ] - cluster_center)^2)
      wcss <- wcss + distance_squared
    }
    return(wcss)
  }
  tryCatch(
    {
      if (k == 1) {
        withinss_record <- c()
        clusters_record <- matrix(nrow = 50, ncol = dim(X)[1])
        for (iter in 1:50) {
          clusters <- rep(1, times = dim(X)[1])
          tkmeans_center <- matrix(apply(X, 2, mean), nrow = 1)
          withinss <- calculate_wcss(X, tkmeans_center, clusters)
          withinss_record <- c(withinss_record, withinss)
          clusters_record[iter, ] <- clusters
        }
        clusters <- clusters_record[which.min(withinss_record), ]
        withinss <- withinss_record[which.min(withinss_record)]
      } else {
        withinss_record <- c()
        clusters_record <- matrix(nrow = 50, ncol = dim(X)[1])
        for (iter in 1:50) {
          init_points <- flexclust::kcca(X, k, flexclust::kccaFamily("kmeans"), control = list(initcent = "kmeanspp"))@centers
          tkmeans_center <- trimcluster::trimkmeans(X, k, trim = 0.1, runs = 1, points = init_points)$means
          clusters <- c(lowmemtkmeans::nearest_cluster(X, tkmeans_center), recursive = TRUE)
          withinss <- calculate_wcss(X, tkmeans_center, clusters)
          withinss_record <- c(withinss_record, withinss)
          clusters_record[iter, ] <- clusters
        }
        clusters <- clusters_record[which.min(withinss_record), ]
        withinss <- withinss_record[which.min(withinss_record)]
      }
    },
    error = function(e) {
      k <- k + 1
      if (k == 1) {
        withinss_record <- c()
        clusters_record <- matrix(nrow = 50, ncol = dim(X)[1])
        for (iter in 1:50) {
          clusters <- rep(1, times = dim(X)[1])
          tkmeans_center <- matrix(apply(X, 2, mean), nrow = 1)
          withinss <- calculate_wcss(X, tkmeans_center, clusters)
          withinss_record <- c(withinss_record, withinss)
          clusters_record[iter, ] <- clusters
        }
        clusters <- clusters_record[which.min(withinss_record), ]
        withinss <- withinss_record[which.min(withinss_record)]
      } else {
        withinss_record <- c()
        clusters_record <- matrix(nrow = 50, ncol = dim(X)[1])
        for (iter in 1:50) {
          init_points <- flexclust::kcca(X, k, flexclust::kccaFamily("kmeans"), control = list(initcent = "kmeanspp"))@centers
          tkmeans_center <- trimcluster::trimkmeans(X, k, trim = 0.1, runs = 1, points = init_points)$means
          clusters <- c(lowmemtkmeans::nearest_cluster(X, tkmeans_center), recursive = TRUE)
          withinss <- calculate_wcss(X, tkmeans_center, clusters)
          withinss_record <- c(withinss_record, withinss)
          clusters_record[iter, ] <- clusters
        }
        clusters <- clusters_record[which.min(withinss_record), ]
        withinss <- withinss_record[which.min(withinss_record)]
      }
    }
  )

  list(cluster = clusters, withinss = withinss)
}

#' @title Community detection and community number selection.
#' @description A function for community detection and community number selection by performing (trimmed) K-means algorithm on fPCA scores.
#' @param X an n*d matrix. fPCA scores obtained from fPCA or fPCA_angle.
#' @param max_k scalar. the maximal number of communities to choose from.
#' @param cluster_method (optional) vector.string. if cluster_method = "k_means", use K-means method for clustering; if cluster_method = "tk_means", use trimmed K-means method for clustering;
#'     defulat "tk_means".
#' @param K_selection (optional) vector.string. if K_selection = "Silhouette", use Silhouette coefficient as the criterion; if K_selection = "Pham", use the method proposed by Pham, Dimov and Nguyen for clustering. default "Silhouette".
#' @param weight (optional) vector. the weight for each data point when clustering.
#' @param parallel (optional) logical. whether to parallelize when selecting the number of communities.
#'
#' @return community membership results with the optimal number of communities.
#' @export
#'
#' @examples
#' data("SBM_example_fASE0.rda")
#' result_SBM0 <- SBM_example_fASE0[[1]]
#' embedding_SBM0 <- fda::eval.fd(1:50, result_SBM0[[1]])
#' PCA_score <- fPCA(embedding_SBM0, 1:50, 2)
#' clus_res <- cluster_no_selection(PCA_score[[1]], 4)
cluster_no_selection <- function(X, max_k, cluster_method = "tk_means", K_selection = "Silhouette", weight = NULL, parallel = TRUE) {
  # X: FPC for K-means or trimmed K-means; fdata for FADP clustering
  set.seed(1)
  if (cluster_method == "k_means" && K_selection == "Silhouette") {
    if (!is.null(weight)) {
      X <- t(weight * t(X))
    }
    result <- factoextra::fviz_nbclust(X, cluster_kmeans, method = "silhouette", k.max = max_k)
    optim_k <- result$data[which.max(result$data$y), "clusters"]
    optim_k <- as.numeric(optim_k)
    optim_result <- list(optim_k, cluster_kmeans(X, optim_k), result$data$y)
  } else if (cluster_method == "k_means" && K_selection == "Pham") {
    if (!is.null(weight)) {
      X <- t(weight * t(X))
    }

    if (parallel) {
      numCores = max(parallel::detectCores()-2,25)
      cl <- parallel::makeCluster(numCores)
      doSNOW::registerDoSNOW(cl)
    }
    result <- kselection::kselection(X, fun_cluster = cluster_kmeans, max_centers = max_k, parallel = parallel, k_threshold = 1)
    optim_k <- result$k
    optim_result <- list(optim_k, cluster_tkmeans(X, optim_k), result$f_k)
    if (parallel) {
      parallel::stopCluster(cl)
    }
  } else if (cluster_method == "tk_means" && K_selection == "Silhouette") {
    if (!is.null(weight)) {
      X <- t(weight * t(X))
    }
    result <- factoextra::fviz_nbclust(X, cluster_tkmeans, method = "silhouette", k.max = max_k)
    optim_k <- result$data[which.max(result$data$y), "clusters"]
    optim_k <- as.numeric(optim_k)
    optim_result <- list(optim_k, cluster_tkmeans(X, optim_k), result$data$y)
  } else if (cluster_method == "tk_means" && K_selection == "Pham") {
    if (!is.null(weight)) {
      X <- t(weight * t(X))
    }
    if (parallel) {
      numCores = max(parallel::detectCores()-2,25)
      cl <- parallel::makeCluster(numCores)
      doSNOW::registerDoSNOW(cl)
    }
    result <- kselection::kselection(X, fun_cluster = cluster_tkmeans, max_centers = max_k, parallel = parallel, k_threshold = 1)
    optim_k <- result$k
    optim_result <- list(optim_k, cluster_tkmeans(X, optim_k), result$f_k)
    if (parallel) {
      parallel::stopCluster(cl)
    }
  }

  optim_result
}

#' @title Visualize adjacency matrices.
#' @description Visualize dynamic networks by adjacency matrices.
#' @param adjacency_matrix matrix. the adjacency matrix of the dynamic network at a specific time point.
#' @param label vector. the community label of nodes.
#' @param permuted logical. If TRUE, the order of nodes in the matrix is permuted according to the label of nodes, and nodes with the same label will be clustered into one group.
#' @param guide_break vector. If NOT NULL, this vector will be used as the aesthetic for shade of colors.
#' @param n the number of total nodes.
#'
#' @return a ggplot2 object.
adjacency_visualize <- function(adjacency_matrix, label = NULL, permuted = FALSE, guide_break = NULL, n = NULL) {
  reperm <- function(accumulated_matrix, index, s, k) {
    if (s > 1 && k > 1) {
      if (s == k) {
        acc1 <- accumulated_matrix[index == s, index == s]
        acc2 <- reperm(accumulated_matrix, index, s, s - 1)
        acc3 <- reperm(accumulated_matrix, index, s - 1, s)
        temp <- cbind(acc1, acc2)
        temp2 <- cbind(acc3, reperm(accumulated_matrix, index, s - 1, s - 1))
        accumulated_matrix_new <- rbind(temp, temp2)
      } else {
        if (s < k) {
          acc <- reperm(accumulated_matrix, index, s - 1, k)
          if (length(which(index == k)) == 1) {
            accumulated_matrix_new <- c(accumulated_matrix[index == s, index == k], acc)
          } else {
            accumulated_matrix_new <- rbind(accumulated_matrix[index == s, index == k], acc)
          }
        } else {
          acc <- reperm(accumulated_matrix, index, s, k - 1)
          if (length(which(index == s)) == 1) {
            accumulated_matrix_new <- cbind(t(accumulated_matrix[index == s, index == k]), acc)
          } else {
            accumulated_matrix_new <- cbind(accumulated_matrix[index == s, index == k], acc)
          }
        }
      }
    } else {
      acc1 <- accumulated_matrix[index == s, index == k]
      if (is.null(dim(acc1)) && length(which(index == s)) == 1) {
        acc1 <- t(acc1)
      }
      accumulated_matrix_new <- acc1
    }
    accumulated_matrix_new
  }


  if (!is.null(label)) {
    k <- range(label)[2]
    t <- rev(table(label))
    t <- cumsum(t)
    if (permuted) adjacency_matrix <- reperm(adjacency_matrix, label, k, k)
  }
  adjacency_matrix[adjacency_matrix == 0] <- NA
  colnames(adjacency_matrix) <- NULL
  permute <- reshape2::melt(adjacency_matrix)
  node_num <- dim(adjacency_matrix)[1]
  permute$Var1 <- node_num - permute$Var1
  permute <- permute[!is.na(permute$value), ]
  if (length(setdiff(1:n, c(permute[, 1:2], recursive = TRUE))) != 0) {
    h <- setdiff(1:n, c(permute[, 1:2], recursive = TRUE))
    d <- rbind(cbind(h, n, NA), cbind(h, 1, NA), cbind(n, h, NA), cbind(1, h, NA))
    colnames(d) <- c("Var1", "Var2", "value")
    permute <- rbind(permute, d)
  }

  p1 <- ggplot(permute, aes(x = Var1, y = Var2, fill = value))
  if (is.null(guide_break)) {
    p2 <- p1 + geom_tile() + scale_x_discrete(name = "Nodes") + scale_y_discrete(name = "Nodes") + scale_fill_gradient2(low = "midnightblue", mid = "cornflowerblue", high = "white", na.value = "transparent")
  } else {
    if (max(guide_break) <= 3) {
      p2 <- p1 + geom_tile() + scale_x_discrete(name = "Nodes") + scale_y_discrete(name = "Nodes") +
        scale_fill_gradientn(colours = c("cornflowerblue", "royalblue", "midnightblue"), breaks = guide_break, labels = guide_break, values = scales::rescale(c(1, max(1 + guide_break) / 2, max(guide_break))), na.value = "transparent")
    } else {
      p2 <- p1 + geom_tile() + scale_x_discrete(name = "Nodes") + scale_y_discrete(name = "Nodes") +
        scale_fill_gradientn(colours = c("cornflowerblue", "royalblue", "midnightblue"), breaks = guide_break, labels = guide_break, values = scales::rescale(c(1, 3, max(guide_break))), na.value = "transparent")
    }
  }
  if (!is.null(label)) {
    p3 <- p2 + geom_segment(x = node_num - t[1] - 0.5, xend = node_num - t[1] - 0.5, y = 0, yend = t[2], color = "pink") + geom_segment(x = node_num - t[2] - 0.5, xend = node_num - 0.5, y = t[1], yend = t[1], color = "pink")
    if (k > 2) {
      for (kk in 2:(k - 1)) {
        p3 <- p3 + geom_segment(x = node_num - t[kk] - 0.5, xend = node_num - t[kk] - 0.5, y = t[kk - 1], yend = t[kk + 1], color = "pink") + geom_segment(x = node_num - t[kk + 1] - 0.5, xend = node_num - t[kk - 1] - 0.5, y = t[kk], yend = t[kk], color = "pink")
      }
    }
    p3 <- p3 + theme(legend.position = "left")
  } else {
    p3 <- p2
  }
  p3
}

#' @title Sort latent functions into a dataframe
#' @description Transform a T*n*d latent function value array into a dataframe for the sake of plotting, without radius correction.
#' @param X a T*n*d latent function value array.
#' @param timestamp_vec T vector. time points corresponding to the first dimension of X.
#'
#' @return a dataframe containing the latent function value for each node at each time point.
#' @export
#'
#' @examples
#' data("SBM_example_fASE0.rda")
#' result_SBM0 <- SBM_example_fASE0[[1]]
#' embedding_SBM0 <- fda::eval.fd(1:50, result_SBM0[[1]])
#' node_embedding_frame0 <- function_sort(embedding_SBM0, 1:50)
function_sort <- function(X, timestamp_vec) {
  node_embedding_frame <- c()
  node_embedding_frame_name <- c()

  if (is.matrix(X)) {
    n_nodes <- dim(X)[2]
    m_layers <- dim(X)[1]

    if (length(timestamp_vec) != m_layers) {
      cat("Parameters do not fit with each other!\n")
      break
    } else {
      for (i in 1:n_nodes) {
        node_embedding_frame <- rbind(node_embedding_frame, cbind(X[, i], t = timestamp_vec, i = i))
      }
      colnames(node_embedding_frame)[1] <- "r"
      node_embedding_frame <- as.data.frame(node_embedding_frame)

      node_embedding_frame
    }
  } else {
    n_nodes <- dim(X)[2]
    m_layers <- dim(X)[1]
    embedding_dim <- dim(X)[3]

    if (length(timestamp_vec) != m_layers) {
      cat("Parameters do not fit with each other!\n")
      break
    } else {
      for (i in 1:n_nodes) {
        node_embedding_frame <- rbind(node_embedding_frame, cbind(X[, i, ], t = timestamp_vec, i = i))
      }
      for (j in 1:embedding_dim) {
        node_embedding_frame_name <- c(node_embedding_frame_name, paste0("x", j))
      }

      colnames(node_embedding_frame)[1:embedding_dim] <- node_embedding_frame_name
      node_embedding_frame <- as.data.frame(node_embedding_frame)

      node_embedding_frame
    }
  }
}

#' @title Sort degree-corrected latent functions into a dataframe
#' @description Transform a T*n*d latent function value array into a dataframe for the sake of plotting, with radius correction.
#' @param X a T*n*d latent function value array.
#' @param timestamp_vec T vector. time points corresponding to the first dimension of X.
#' @param active_points list. a list of time points for each node. if NULL, regard all nodes active at all time points. default NULL.
#'
#' @return a dataframe containing the latent function value for each node at each time point.
#' @export
#'
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
#' ER_varynodes <- function(n, p, t) {
#'   adjacency_matrix <- matrix(0, nrow = n, ncol = n)
#'
#'   active_nodes <- unique(c((floor((t - 1) / 10) * 10 + 1):((floor((t - 1) / 10) + 1) * 10), 80:100))
#'   for (i in active_nodes) {
#'     for (j in active_nodes) {
#'       if (i > j) {
#'         adjacency_matrix[i, j] <- sample(c(0, 1), size = 1, prob = c(1 - p, p))
#'       }
#'     }
#'     adjacency_matrix[i, i] <- 0.5
#'   }
#'   adjacency_matrix + t(adjacency_matrix)
#' }

#' set.seed(10)
#' n_nodes = 100
#' T = 100
#' p = 0.1
#' C = c(rep(2,times=79),rep(1,times=21))
#' dynamic_adjacency = c()
#' dynamic_network_adjacency = array(0, dim = c(n_nodes, n_nodes,T))
#' for(t in 1:T){
#'   adj = ER_varynodes(n_nodes,p,t)
#'   dynamic_adjacency = c(dynamic_adjacency,list(list(adj)))
#'   dynamic_network_adjacency[,,t] = adj
#' }
#' data("varySBM_example_fASE0.rda")
#' result_varySBM0 = varySBM_example_fASE0[[1]]
#' embedding_DCBM0 = fda::eval.fd(1:100, result_varySBM0[[1]])
#' active_points = lapply(1:100, active_calculation, dynamic_network_adjacency, timestamp_vec = 1:100)
#' node_embedding_frame0 = function_sort_angle(embedding_DCBM0, 1:100,active_points=active_points)
#' node_embedding_frame0 = cbind(node_embedding_frame0, label = factor(rep(C, each = T)))
#' func_plot_x1 = ggplot(node_embedding_frame0,aes_string(x="t",y=paste0("x",1),group="i",color="label"))+geom_line(aes(alpha=label,size=label))+scale_color_manual(values=c("#FA7F6F","#82B0D2"))+scale_size_manual(values= c(0.8,0.8))+scale_alpha_manual(values= c(0.8,0.8))+guides(color=guide_legend(override.aes = list(size=2)))+theme_bw()
function_sort_angle <- function(X, timestamp_vec, active_points = NULL) {
  node_embedding_frame <- c()
  node_embedding_frame_name <- c()
  subset <- function(i, ls, x) {
    if (length(ls[[i]]) != 0) {
      h <- matrix(0, nrow = dim(x)[1], ncol = dim(x)[3])
      h[which((1:length(timestamp_vec)) %in% ls[[i]]), ] <- x[which((1:length(timestamp_vec)) %in% ls[[i]]), i, ]
    } else {
      h <- matrix(0, nrow = dim(x)[1], ncol = dim(x)[3])
    }
    h
  }
  subset_mat <- function(i, ls, x) {
    if (length(ls[[i]]) != 0) {
      h <- matrix(0, nrow = dim(x)[1], ncol = 1)
      h[which((1:length(timestamp_vec)) %in% ls[[i]]), ] <- x[which((1:length(timestamp_vec)) %in% ls[[i]]), i]
    } else {
      h <- matrix(0, nrow = dim(x)[1], ncol = 1)
    }
    h
  }

  if (!is.null(active_points)) {
    active_points <- lapply(active_points, function(x) {
      match(x, timestamp_vec)
    })
  }

  if (is.matrix(X)) {
    n_nodes <- dim(X)[2]
    m_layers <- dim(X)[1]

    embedding_fASE2 <- X
    if (!is.null(active_points)) {
      embedding_fASE2 <- sapply(1:dim(X)[2], subset_mat, active_points, embedding_fASE2, simplify = "array")
    }
    X <- embedding_fASE2

    if (!is.null(active_points)) {
      active_number <- sapply(1:dim(X)[2], function(x) {
        length(active_points[[x]])
      })
      active_number[active_number == 0] <- 1
      X <- t(t(X) / sqrt(active_number))
    }


    if (length(timestamp_vec) != m_layers) {
      cat("Parameters do not fit with each other!\n")
      break
    } else {
      for (i in 1:n_nodes) {
        node_embedding_frame <- rbind(node_embedding_frame, cbind(X[, i], t = timestamp_vec, i = i))
      }
      colnames(node_embedding_frame)[1] <- "r"
      node_embedding_frame <- as.data.frame(node_embedding_frame)

      node_embedding_frame
    }
  } else {
    n_nodes <- dim(X)[2]
    m_layers <- dim(X)[1]
    embedding_dim <- dim(X)[3]

    embedding_fASE2 <- X
    if (!is.null(active_points)) {
      embedding_fASE2 <- sapply(1:dim(X)[2], subset, active_points, embedding_fASE2, simplify = "array")
      embedding_fASE2 <- aperm(embedding_fASE2, c(1, 3, 2))
    }
    embedding_fASE3 <- apply(embedding_fASE2, c(1, 2), crossprod)
    embedding_fASE3[embedding_fASE3 == 0] <- 1

    X <- embedding_fASE2 / sqrt(array(embedding_fASE3, dim = c(dim(X)[1], dim(X)[2], dim(X)[3])))
    if (!is.null(active_points)) {
      active_number <- sapply(1:dim(X)[2], function(x) {
        length(active_points[[x]])
      })
      active_number[active_number == 0] <- 1
      X <- X / sqrt(aperm(array(active_number, dim = c(dim(X)[2], dim(X)[1], dim(X)[3])), c(2, 1, 3)))
    }

    if (length(timestamp_vec) != m_layers) {
      cat("Parameters do not fit with each other!\n")
      break
    } else {
      for (i in 1:n_nodes) {
        node_embedding_frame <- rbind(node_embedding_frame, cbind(X[, i, ], t = timestamp_vec, i = i))
      }
      for (j in 1:embedding_dim) {
        node_embedding_frame_name <- c(node_embedding_frame_name, paste0("x", j))
      }

      colnames(node_embedding_frame)[1:embedding_dim] <- node_embedding_frame_name
      node_embedding_frame <- as.data.frame(node_embedding_frame)

      node_embedding_frame
    }
  }
}

#' @title draw the cd(community detection)/dccd(degree-corrected community detection) plot
#' @description Perform a (radius-corrected) MVFPCA on latent functions, and then draw the (degree-corrected) community detection plot.
#' @param embedding_fASE an T*n*d array, containing values of the resulting latent functions obtained by fASE function at timestamp_vec.
#' @param timestamp_vec T vector. the observed time points of observation.
#' @param radius_corrected (optional) logical. if radius_corrected = TRUE, then each latent function is divided by its radius function to correct the degree effect. default TRUE.
#' @param active_corrected (optional) logical. if active_corrected = TRUE, then each latent function is divided by the active time of the corresponding node to correct the active duration effect. default FALSE.
#' @param active_points (optional) list whose length is the number of nodes. when active_corrected = TRUE, active_points is required to provide the active time of each node.
#' @param cluster (optional) vector. the community label obtained by some clustering method.
#' @param label (optional) vector. the true community label.
#' @param cluster_color (optional) vector. colors for each community obtained by some clustering method.
#' @param label_color  (optional) vector. colors for each true community.
#' @param point_size (optional) scalar. point size.
#'
#' @return a ggplot2 object of the cd plot.
#' @export
#'
#' @examples
#' data("DCBM_example_fASE0.rda")
#' result_DCBM0 <- DCBM_example_fASE0[[1]]
#' deg <- DCBM_example_fASE0[[2]]
#' embedding_DCBM0 <- fda::eval.fd(1:100, result_DCBM0[[1]])
#' PCA_score <- fPCA_angle(embedding_DCBM0, 1:100, 20)
#' clus_res <- cluster_no_selection(PCA_score[[1]], 10)
#' PCA_score <- fPCA_angle(embedding_DCBM0, 1:100, 2)
#' C <- c(rep(1, times = 25), rep(3, times = 25), rep(5, times = 25), rep(2, times = 25), rep(4, times = 10))
#' cdplot <- cd_plot(embedding_DCBM0, 1:100, cluster = clus_res[[2]]$cluster, label = C)
cd_plot <- function(embedding_fASE, timestamp_vec, radius_corrected = TRUE, active_corrected = FALSE, active_points = NULL, cluster = NULL, label = NULL, cluster_color = NULL, label_color = NULL, point_size = NULL) {
  if (radius_corrected) {
    if (active_corrected) {
      PCA_score <- fPCA_angle(embedding_fASE, timestamp_vec, 2, active_points)
    } else {
      PCA_score <- fPCA_angle(embedding_fASE, timestamp_vec, 2)
    }
    PCA_score <- as.data.frame(PCA_score[[1]])
    names(PCA_score) <- c("theta1", "theta2")
  } else {
    if (active_corrected) {
      PCA_score <- fPCA(embedding_fASE, timestamp_vec, 2, active_points)
    } else {
      PCA_score <- fPCA(embedding_fASE, timestamp_vec, 2)
    }
    PCA_score <- as.data.frame(PCA_score[[1]])
    names(PCA_score) <- c("x1", "x2")
  }
  if (!is.null(cluster)) {
    PCA_score$cluster <- factor(cluster)
  }
  if (!is.null(label)) {
    PCA_score$label <- factor(label)
  }
  if (is.null(cluster_color)) {
    cluster_color <- c("#8ECFC9", "#FFBE7A", "#FA7F6F", "#82B0D2", "#BEB8DC", "#E7DAD2", "gray")
  }
  if (is.null(label_color)) {
    label_color <- c("#8ECFC9", "#FFBE7A", "#FA7F6F", "#82B0D2", "#BEB8DC", "#E7DAD2", "gray")
  }
  if (is.null(point_size)) {
    point_size <- 3
  }

  if (!is.null(cluster) && !is.null(label)) {
    plot <- ggplot(PCA_score, aes(x = theta1, y = theta2, group = cluster, color = label)) +
      geom_point() +
      ggalt::geom_encircle(aes(group = cluster, fill = cluster), expand = 0.03, spread = 0.5, s_shape = 1, size = point_size, linetype = 1, alpha = 0.2) +
      scale_color_manual(values = label_color) +
      scale_fill_manual(values = cluster_color) +
      guides(fill = "none") +
      theme_bw()
  } else if (is.null(cluster) && !is.null(label)) {
    plot <- ggplot(PCA_score, aes(x = theta1, y = theta2, color = label)) +
      geom_point() +
      scale_color_manual(values = label_color) +
      guides(fill = "none") +
      theme_bw()
  } else if (!is.null(cluster) && is.null(label)) {
    plot <- ggplot(PCA_score, aes(x = theta1, y = theta2, group = cluster)) +
      geom_point() +
      ggalt::geom_encircle(aes(group = cluster, fill = cluster), expand = 0.03, spread = 0.5, s_shape = 1, size = point_size, linetype = 1, alpha = 0.2) +
      scale_fill_manual(values = cluster_color) +
      guides(fill = "none") +
      theme_bw()
  } else {
    plot <- ggplot(PCA_score, aes(x = theta1, y = theta2)) +
      geom_point() +
      guides(fill = "none") +
      theme_bw()
  }

  plot
}


#' @title draw the dd (degree-depth) plot
#' @description Obtain the statistical depth for each node by calculating the cosine distances between latent functions, and draw the degree-depth plot.
#' @param embedding_fASE a T*n*d array, containing values of the resulting latent functions obtained by fASE function at timestamp_vec.
#' @param timestamp_vec T vector. the observed time points of observation.
#' @param degree_vec n vector. the total degree of each node in the network.
#' @param cluster  (optional) vector. the community label obtained by some clustering method.
#' @param label (optional) vector. the true community label.
#' @param cluster_shape (optional) vector. shapes for each community obtained by some clustering method.
#' @param label_color  (optional) vector. colors for each true community.
#' @param point_size (optional) scalar. point size.
#'
#' @return  a ggplot2 object of the dd plot.
#' @export
#'
#' @examples
#' data("DCBM_example_fASE0.rda")
#' result_DCBM0 <- DCBM_example_fASE0[[1]]
#' deg <- DCBM_example_fASE0[[2]]
#' embedding_DCBM0 <- fda::eval.fd(1:100, result_DCBM0[[1]])
#' PCA_score <- fPCA_angle(embedding_DCBM0, 1:100, 20)
#' clus_res <- cluster_no_selection(PCA_score[[1]], 10)
#' C <- c(rep(1, times = 25), rep(3, times = 25), rep(5, times = 25), rep(2, times = 25), rep(4, times = 10))
#' ddplot <- dd_plot(embedding_DCBM0, 1:100, deg, cluster = clus_res[[2]]$cluster, label = C)
dd_plot <- function(embedding_fASE, timestamp_vec, degree_vec, cluster = NULL, label = NULL, cluster_shape = NULL, label_color = NULL, point_size = NULL) {
  sim <- cosine_similarity(embedding_fASE, timestamp_vec)
  depth_stats <- metric_depth(1 - sim)
  depth_reindex <- order(depth_stats, decreasing = TRUE)

  dd_frame <- data.frame(depth = depth_stats)
  dd_frame[depth_reindex, "depth_order"] <- as.numeric(1:length(depth_stats))
  dd_frame$degree <- degree_vec

  if (!is.null(cluster)) {
    dd_frame$cluster <- factor(cluster)
  }
  if (!is.null(label)) {
    dd_frame$label <- factor(label)
  }
  if (is.null(cluster_shape)) {
    cluster_shape <- c(15, 17, 16, 18, 1:14)
  }
  if (is.null(label_color)) {
    label_color <- c("#8ECFC9", "#FFBE7A", "#FA7F6F", "#82B0D2", "#BEB8DC", "#E7DAD2", "gray")
  }
  if (is.null(point_size)) {
    point_size <- 3
  }

  if (!is.null(cluster) && !is.null(label)) {
    ddplot <- ggplot(dd_frame, aes(x = depth, y = degree, color = label)) +
      geom_point(aes(shape = cluster), size = point_size) +
      scale_y_log10() +
      annotation_custom(
        grob = grid::linesGrob(y = unit(c(1, 0), "npc")), xmin = -Inf, ymin = Inf, xmax = Inf, ymax = -Inf
      ) +
      scale_shape_manual(values = cluster_shape) +
      scale_color_manual(values = label_color) +
      guides(color = "none", fill = "none", shape = "none") +
      theme_bw() +
      theme(text = element_text(size = 28), axis.title = element_text(size = 24), legend.text = element_text(size = 24))
  } else if (is.null(cluster) && !is.null(label)) {
    ddplot <- ggplot(dd_frame, aes(x = depth, y = degree, color = label)) +
      geom_point(aes(shape = cluster), size = point_size) +
      scale_y_log10() +
      annotation_custom(
        grob = grid::linesGrob(y = unit(c(1, 0), "npc")), xmin = -Inf, ymin = Inf, xmax = Inf, ymax = -Inf
      ) +
      scale_color_manual(values = label_color) +
      guides(color = "none", fill = "none", shape = "none") +
      theme_bw() +
      theme(text = element_text(size = 28), axis.title = element_text(size = 24), legend.text = element_text(size = 24))
  } else if (!is.null(cluster) && is.null(label)) {
    ddplot <- ggplot(dd_frame, aes(x = depth, y = degree)) +
      geom_point(aes(shape = cluster), size = point_size) +
      scale_y_log10() +
      annotation_custom(
        grob = grid::linesGrob(y = unit(c(1, 0), "npc")), xmin = -Inf, ymin = Inf, xmax = Inf, ymax = -Inf
      ) +
      scale_shape_manual(values = cluster_shape) +
      guides(color = "none", fill = "none", shape = "none") +
      theme_bw() +
      theme(text = element_text(size = 28), axis.title = element_text(size = 24), legend.text = element_text(size = 24))
  } else {
    ddplot <- ggplot(dd_frame, aes(x = depth, y = degree)) +
      geom_point(size = point_size) +
      scale_y_log10() +
      annotation_custom(
        grob = grid::linesGrob(y = unit(c(1, 0), "npc")), xmin = -Inf, ymin = Inf, xmax = Inf, ymax = -Inf
      ) +
      guides(color = "none", fill = "none", shape = "none") +
      theme_bw() +
      theme(text = element_text(size = 28), axis.title = element_text(size = 24), legend.text = element_text(size = 24))
  }

  ddplot
}
