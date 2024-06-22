test_that("accelerated fASE works", {
  ER_generation <- function(n, p) {
    adjacency_matrix <- matrix(0, nrow = n, ncol = n)
    for (i in 2:n) {
      for (j in 1:(i - 1)) {
        adjacency_matrix[i, j] <- sample(c(0, 1), size = 1, prob = c(1 - p, p))
      }
      adjacency_matrix[i, i] <- 0.5
    }
    adjacency_matrix + t(adjacency_matrix)
  }
  set.seed(1)
  n_nodes <- 10
  T <- 10
  p <- 0.2
  dynamic_adjacency <- c()
  dynamic_network_adjacency <- array(0, dim = c(n_nodes, n_nodes, T))
  for (t in 1:T) {
    adj <- ER_generation(n_nodes, p)
    dynamic_adjacency <- c(dynamic_adjacency, list(list(adj)))
    dynamic_network_adjacency[, , t] <- adj
  }
  result_ER1 <- fASE(dynamic_network_adjacency, 2, 4, epsilon = 1, batch_size = NULL)
  expect_identical(result_ER1, result_ER2)
})

test_that("cd plot works", {
  result_DCBM0 <- DCBM_example_fASE0[[1]]
  deg <- DCBM_example_fASE0[[2]]
  embedding_DCBM0 <- fda::eval.fd(1:100, result_DCBM0[[1]])
  PCA_score <- fPCA_angle(embedding_DCBM0, 1:100, 20)
  clus_res <- cluster_no_selection(PCA_score[[1]], 10)
  PCA_score <- fPCA_angle(embedding_DCBM0, 1:100, 2)
  C <- c(rep(1, times = 25), rep(3, times = 25), rep(5, times = 25), rep(2, times = 25), rep(4, times = 10))
  cdplot <- cd_plot(embedding_DCBM0, 1:100, cluster = clus_res[[2]]$cluster, label = C)
  equivalent_ggplot <- function(x, y) {
    # ggplot_table triggers a blank plot that can't be silenced so we divert it
    # not sure if pdf() is the most efficient
    pdf(tempfile(fileext = ".pdf"))
    x_tbl <- suppressWarnings(ggplot2::ggplot_gtable(ggplot2::ggplot_build(x)))
    y_tbl <- suppressWarnings(ggplot2::ggplot_gtable(ggplot2::ggplot_build(y)))
    dev.off()
    # we could probably do a better index equivalency check than just scrubbing
    # them off, but I haven't figured out how it works
    x_unlisted <- gsub("\\d+", "XXX", unlist(x_tbl))
    y_unlisted <- gsub("\\d+", "XXX", unlist(y_tbl))
    names(x_unlisted) <- gsub("\\d+", "XXX", names(x_tbl))
    names(y_unlisted) <- gsub("\\d+", "XXX", names(y_tbl))
    expect_identical(x_unlisted, y_unlisted)
  }
  equivalent_ggplot(cdplot, cdres)
})

test_that("dd plot works", {
  result_DCBM0 <- DCBM_example_fASE0[[1]]
  deg <- DCBM_example_fASE0[[2]]
  embedding_DCBM0 <- fda::eval.fd(1:100, result_DCBM0[[1]])
  PCA_score <- fPCA_angle(embedding_DCBM0, 1:100, 20)
  clus_res <- cluster_no_selection(PCA_score[[1]], 10)
  C <- c(rep(1, times = 25), rep(3, times = 25), rep(5, times = 25), rep(2, times = 25), rep(4, times = 10))
  ddplot <- dd_plot(embedding_DCBM0, 1:100, deg, cluster = clus_res[[2]]$cluster, label = C)
  equivalent_ggplot <- function(x, y) {
    # ggplot_table triggers a blank plot that can't be silenced so we divert it
    # not sure if pdf() is the most efficient
    pdf(tempfile(fileext = ".pdf"))
    x_tbl <- suppressWarnings(ggplot2::ggplot_gtable(ggplot2::ggplot_build(x)))
    y_tbl <- suppressWarnings(ggplot2::ggplot_gtable(ggplot2::ggplot_build(y)))
    dev.off()
    # we could probably do a better index equivalency check than just scrubbing
    # them off, but I haven't figured out how it works
    x_unlisted <- gsub("\\d+", "XXX", unlist(x_tbl))
    y_unlisted <- gsub("\\d+", "XXX", unlist(y_tbl))
    names(x_unlisted) <- gsub("\\d+", "XXX", names(x_tbl))
    names(y_unlisted) <- gsub("\\d+", "XXX", names(y_tbl))
    expect_identical(x_unlisted, y_unlisted)
  }
  equivalent_ggplot(ddplot, ddres)
})


