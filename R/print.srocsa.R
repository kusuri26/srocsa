#' @title
#' Print method for 'srocsa' objects
#'
#' @description This function prints the main results:
#' convergence of the optimizer, SAUC estimates, and parameter estimates with standard errors (SE)
#' for each marginal selection probability \eqn{p}.
#' For the time-dependent results, the results are printed at each time point.
#'
#' @param x
#' A fitted object returned by the function \code{srocsa.t}, \code{tsroc}, or \code{tsrocsa.lrt}.
#'
#' @param digits
#' The digits of the results
#'
#' @param ... Other arguments in function \code{\link{print.default}}
#'
#' @details
#' \describe{
#' \item{Optimizer}{1 indicates that \code{nlminb} converged (\code{convergence = 0}).}
#' \item{SAUC}{The SE is calculated by the delta method.
#' The confidence interval is calculated on the logit scale, so it is not symmetric around the SAUC.}
#' \item{tau_se, tau_sp}{Between-study standard deviations (not variances).}
#' \item{beta, alpha}{"-" indicates not applicable (\eqn{p = 1}, or the model without selection).
#' The SE of alpha is not available.}
#' \item{Results from \code{tsroc}}{Shown as \eqn{p = 1}, because \code{tsroc} fits the model without selection.}
#' }
#'
#' @examples
#' ## load data
#' data(IVD)
#'
#' ## print the results
#' fit1 <- srocsa.t(TP = TP, FN = FN, TN = TN, FP = FP, data = IVD, sa.p = c(1, 0.8, 0.6))
#' print(fit1)
#'
#' fit2 <- srocsa.t(TP = TP, FN = FN, TN = TN, FP = FP, data = IVD, sa.p = c(1, 0.8, 0.6), c11 = 0.5)
#' print(fit2, digits = 4)
#'
#' @seealso
#' \code{\link{srocsa.t}},
#' \code{\link{tsroc}},
#' \code{\link{tsrocsa.lrt}},
#' \code{\link{summary.srocsa}},
#' \code{\link{plot.srocsa}}
#'
#' @method print srocsa
#' @return Print the main results of a srocsa object
#' @export

print.srocsa <- function(x, digits = 3, ...){

  blocks <- .srocsa.blocks(x)

  for (b in blocks) .srocsa.print.block(b, digits = digits)

  invisible(x)

}
