#' @docType data
#' @name SBM_example_fASE0
#'
#' @title An illustrative example where snapshots are independent SBMs.
#' @description  A dataset containing the results of accelerated fASE, where snapshots are independent SBMs.
#'
#' @format A list containing a data frame and a list:
#' \describe{
#' \item{result_SBM0}{a fd object containing the latent functions of accelerated fASE, including dunction fd, coefficient array, basis fd, basis splines}
#' \item{tuning_parameter0}{The hyperparameter tuning results.}
#' }
#' @export
"SBM_example_fASE0"

#' @docType data
#' @name DCBM_example_fASE0
#'
#' @title An illustrative example where snapshots are independent DCBMs, overlaid with a core-periphery structure.
#' @description  A dataset containing the results of accelerated fASE, where snapshots are independent SBMs.
#'
#' @format A list containing a data frame, a vector and a list
#' \describe{
#'   \item{result_DCBM0}{a fd object containing the latent functions of accelerated fASE, including dunction fd, coefficient array, basis fd, basis splines}
#'   \item{deg}{The degree of nodes in the dynamic network.}
#'   \item{tuning_parameter0}{The hyperparameter tuning results.}
#' }
#' @export
"DCBM_example_fASE0"

#' @docType data
#' @name varySBM_example_fASE0
#'
#' @title An illustrative example where snapshots are independent DCBMs, overlaid with a core-periphery structure.
#' @description  A dataset containing the results of accelerated fASE, where snapshots are independent SBMs.
#'
#' @format A list containing a data frame and a vector:
#' \describe{
#'   \item{result_varySBM0}{a fd object containing the latent functions of accelerated fASE, including dunction fd, coefficient array, basis fd, basis splines}
#'   \item{C}{The community label of nodes in the dynamic network.}
#' }
#' @export
"varySBM_example_fASE0"
