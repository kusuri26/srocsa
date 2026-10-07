#' @title
#' Plot method for 'srocsa' objects
#'
#' @description
#' This function plots the SROC curves and the summary points for each marginal selection probability \eqn{p}.
#' For the time-dependent results, one plot is created at each time point, with the title \code{time = xx}.
#'
#' @param x
#' A fitted object returned by the function \code{srocsa.t}, \code{tsroc}, or \code{tsrocsa.lrt}.
#'
#' @param time
#' Time points to plot (only for the results from \code{tsroc} or \code{tsrocsa.lrt}).
#' If \code{time=NULL} (default), all the time points are plotted.
#'
#' @param sroc.cols Colors of the SROC curves, one for each \eqn{p}.
#' Default uses grey's colors.
#'
#' @param sroc.lty Line type of SROC curves.
#'
#' @param sroc.lwd Line width of SROC curves.
#'
#' @param add.spoint Whether to add the summary points on the SROC curves.
#'
#' @param sp.pch Type of the summary points.
#'
#' @param sp.cex Size of the summary points.
#'
#' @param legend.pos Position of the legend. See \code{\link[graphics]{legend}}.
#'
#' @param main Title of the plot.
#' If \code{main=NULL} (default), \code{time = xx} is used for the time-dependent results,
#' and no title for the results from \code{srocsa.t}.
#'
#' @param xlab Label of x-axis.
#'
#' @param ylab Label of y-axis.
#'
#' @param ... Other arguments in function \code{\link{plot.default}} or function \code{\link{curve}}
#'
#' @importFrom grDevices gray.colors
#' @importFrom graphics legend title
#'
#' @examples
#' ## load data
#' data(IVD)
#'
#' fit1 <- srocsa.t(TP = TP, FN = FN, TN = TN, FP = FP, data = IVD, sa.p = c(1, 0.8, 0.6))
#' plot(fit1)
#'
#' ## set the colors
#' plot(fit1, sroc.cols = c("black", "blue", "red"), legend.pos = "bottomright")
#'
#' @seealso
#' \code{\link{plot.sroc}},
#' \code{\link{srocsa.t}},
#' \code{\link{tsroc}},
#' \code{\link{tsrocsa.lrt}}
#'
#' @method plot srocsa
#' @return Plot the SROC curves from a srocsa object
#' @export

plot.srocsa <- function(
  x,
  time = NULL,
  sroc.cols = NULL,
  sroc.lty = 1,
  sroc.lwd = 1,
  add.spoint = TRUE,
  sp.pch = 19,
  sp.cex = 1,
  legend.pos = "bottomright",
  main = NULL,
  xlab = "FPR",
  ylab = "TPR",
  ...
){

  blocks <- .srocsa.blocks(x)

  ## select time points ----
  if (!is.null(time)) {

    times <- vapply(blocks, function(b) if (is.null(b$time)) NA_character_ else b$time, character(1))

    if (all(is.na(times))) stop("`time` can be used only for the results from tsroc or tsrocsa.lrt.", call. = FALSE)

    keep <- times %in% as.character(time)
    if (!any(keep)) stop("No results at time = ", paste(time, collapse = ", "), call. = FALSE)

    blocks <- blocks[keep]
  }

  if (!is.null(main)) main <- rep_len(main, length(blocks))

  for (i in seq_along(blocks)) {

    b <- blocks[[i]]
    k <- length(b$p)

    ## one color for each p (the same p has the same color at all time points)
    cols <- if (is.null(sroc.cols)) gray.colors(k, gamma = 1, start = 0, end = 0.5) else rep_len(sroc.cols, k)

    ## skip the curves that could not be estimated
    par.mat <- b$est[, c("mu1", "mu2", "tau1", "tau2", "rho"), drop = FALSE]
    ok      <- complete.cases(par.mat)

    if (!any(ok)) {
      warning("No SROC curve to plot", if (!is.null(b$time)) paste0(" at time = ", b$time), call. = FALSE)
      next
    }

    plot.sroc(
      t(par.mat[ok, , drop = FALSE]),
      sroc.cols  = cols[ok],
      sroc.lty   = sroc.lty,
      sroc.lwd   = sroc.lwd,
      add.spoint = add.spoint,
      sp.pch     = sp.pch,
      sp.cex     = sp.cex,
      xlab       = xlab,
      ylab       = ylab,
      ...
    )

    title(main = if (!is.null(main)) main[i] else if (!is.null(b$time)) paste0("time = ", b$time))

    legend(
      legend.pos,
      legend = paste0("p = ", b$p[ok]),
      col    = cols[ok],
      lty    = sroc.lty,
      lwd    = sroc.lwd,
      pch    = if (add.spoint) sp.pch else NA,
      bty    = "n"
    )
  }

  invisible(x)

}
