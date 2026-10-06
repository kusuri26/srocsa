#' @title Summary method for 'srocsa' objects
#'
#' @description This function summarizes the results:
#' the main results (see \code{\link{print.srocsa}}) and all the elements of the fitted models
#' at each time point and each marginal selection probability \eqn{p}.
#'
#' @param object
#' A fitted object returned by the function \code{srocsa.t}, \code{tsroc}, or \code{tsrocsa.lrt}.
#'
#' @param digits
#' The digits of the results
#'
#' @param ... Other arguments in function \code{\link{summary}}
#'
#' @examples
#' ## load data
#' data(IVD)
#'
#' fit1 <- srocsa.t(TP = TP, FN = FN, TN = TN, FP = FP, data = IVD, sa.p = c(1, 0.8, 0.6))
#' summary(fit1)
#'
#' ## tables of the estimates
#' res <- summary(fit1)
#' res$sauc
#' res$coef
#'
#' @seealso
#' \code{\link{srocsa.t}},
#' \code{\link{tsroc}},
#' \code{\link{tsrocsa.lrt}},
#' \code{\link{print.srocsa}},
#' \code{\link{print.summary.srocsa}}
#'
#' @method summary srocsa
#' @return A list that consists of summaries of a srocsa object:
#' call, type of the model, table of SAUC estimates (\code{sauc}),
#' table of parameter estimates and SE (\code{coef}),
#' and the fitted models at each time point and each \eqn{p} (\code{blocks})
#' @export

summary.srocsa <- function(object, digits = 3, ...) {

  blocks <- .srocsa.blocks(object)

  ## stack the estimates of all time points ----
  stack <- function(cols) {
    tab <- do.call(rbind, lapply(blocks, function(b) {
      d <- data.frame(p = b$p, b$est[, cols, drop = FALSE], row.names = NULL)
      if (!is.null(b$time)) d <- cbind(time = as.numeric(b$time), d)
      d
    }))
    rownames(tab) <- NULL
    tab
  }

  res <- list(
    call   = object[["call"]],
    type   = .srocsa.type(object),
    sauc   = stack(c("conv", "sauc", "sauc.se", "sauc.lb", "sauc.ub")),
    coef   = stack(c("mu1", "mu1.se", "mu2", "mu2.se", "tau1", "tau1.se", "tau2", "tau2.se",
                     "rho", "rho.se", "beta", "beta.se", "alpha", "alpha.se")),
    blocks = blocks
  )
  class(res) <- "summary.srocsa"

  print(res, digits = digits)

  invisible(res)

}
