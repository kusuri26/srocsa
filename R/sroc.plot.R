
#' @title Summary ROC curves
#'
#' @description Plot the SROC curves
#'
#' @param object (1) The result from function \code{dtametasa.fc} or \code{dtametasa.rc} (single curve);
#' (2) a matrix of which the rows are \code{c(u1, u2, t1, t2, r)} (multiple curves);
#' (3) a vector of \code{c(u1, u2, t1, t2, r)} (single curve).
#'
#' @param sroc.type Plot Reitsma's SROC curve(\code{"sroc"}) or Rutter's HSROC curve (\code{"hsroc"}),
#'
#' @param sroc.cols Set a vector of colors for single/multiple SROC curves.
#' Default uses grey's colors.
#'
#' @param sroc.lty Line type of SROC curve.
#'
#' @param sroc.lwd Line width of SROC curve.
#'
#' @param add.spoint Whether to add the summary point on the SROC curve.
#'
#' @param sp.pch Type of the summary point.
#'
#' @param sp.cex Size of the summary point.
#'
#' @param plot.ci Whether to plot the confidence intervals of the SROC curve.
#' It is only for the result from function \code{dtametasa.fc} or \code{dtametasa.rc}
#'
#' @param sroc.ci.col Colors of the confidence intervals.
#'
#' @param sroc.ci.lty Line type of the confidence intervals.
#'
#' @param sroc.ci.lwd Line width of the confidence intervals.
#'
#' @param ci.level Significant level of the confidence intervals.
#'
#' @param xlab Label of x-axis.
#'
#' @param ylab Label of y-axis.
#'
#' @param addon Whether to add the SROC curves onto an existed plot.
#' Default is \code{FALSE}, to start from a new plot.
#'
#' @param ... Other augments in function \code{\link{plot.default}} or function \code{\link{curve}}
#'
#' @importFrom grDevices gray.colors
#' @importFrom graphics curve points
#'
#' @return SROC curves
#'
#' @seealso
#' \code{\link[graphics]{plot.default}},
#' \code{\link[graphics]{points}},
#' \code{\link[graphics]{curve}},
#' \code{\link{dtametasa.fc}},
#' \code{\link{dtametasa.rc}}.
#'
#'
#' @examples
#' p.seq <- seq(1, 0.5, -0.1)
#' sa1.seq <- sapply(p.seq, function (p) dtametasa.fc(IVD, p)$par)
#'
#' SROC(sa1.seq, sroc.type = "sroc")
#' SROC(sa1.seq, sroc.type = "hsroc")
#'
#' sa.fit <- dtametasa.fc(IVD, p = 0.7, sauc.type = "sroc")
#' sa.fit
#' SROC(sa.fit, sroc.type = "sroc", pch = 19, plot.ci = TRUE)
#' SROC(sa.fit, sroc.type = "hsroc", pch = 19, plot.ci = TRUE)
#'
#' @export

sroc <- function(
  # object,
  # sroc.type = c("sroc", "hsroc"),
  # mu1, mu2, tau1, tau2, rho,
  par.mat,
  var.matrix = NULL,
  sroc.cols = gray.colors(ncol(as.matrix(par.mat)), gamma = 1, start = 0, end = 0.5),
  sroc.lty = 1,
  sroc.lwd = 1,
  add.spoint = TRUE,
  sp.pch = 19,
  sp.cex = 1,
  plot.ci = FALSE,
  sroc.ci.col = "red",
  sroc.ci.lty = 2,
  sroc.ci.lwd = 1,
  ci.level = 0.95,
  xlab = "FPR",
  ylab = "TPR",
  addon = FALSE,
  ...
){

  par <- as.matrix(par.mat)

  if(nrow(par) < 5) stop("PLEASE CHECK THE INPUT MATRIX")

  # sroc.type <- match.arg(sroc.type)

  if (!addon) plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab = xlab, ylab = ylab, ...)

  for (i in 1:ncol(par)) {

    u1  <- par[1,i]
    u2  <- par[2,i]
    t1  <- par[3,i]
    t2  <- par[4,i]
    r <- par[5, i] 

    f <- function(x) plogis(u1 - (t1*t2*r/(t2^2)) * (qlogis(x) + u2))
    curve(f, 0, 1,
          col = sroc.cols[i],
          add = TRUE,
          lty = sroc.lty,
          lwd = sroc.lwd,
          ...)
  }

  if (add.spoint) {

    sens <- plogis(par[1,])
    spec <- plogis(par[2,])
    points(1-spec, sens, col=sroc.cols, pch = sp.pch, cex = sp.cex)

    }

  if(plot.ci && !is.null(var.matrix)){

    for (i in 1:ncol(par)) {

    u1 <- par[1,]
    u2 <- par[2,]
    t1 <- par[3,]
    t2 <- par[4,]
    r  <- par[5,]


      f.lb <- function(x) plogis( u1 - (t1*t2*r/t2^2) * (qlogis(x) + u2) + qnorm((1-ci.level)/2, lower.tail = TRUE) * suppressWarnings(sqrt(.QIQ.sroc(x, mu1=u1, mu2=u2, tau1=t1, tau2=t2, rho=r, var.matrix[1:5,1:5]))))

      f.ub <- function(x) plogis( u1 - (t1*t2*r/t2^2) * (qlogis(x) + u2) + qnorm((1-ci.level)/2, lower.tail = FALSE) * suppressWarnings(sqrt(.QIQ.sroc(x, mu1=u1, mu2=u2, tau1=t1, tau2=t2, rho=r, var.matrix[1:5,1:5]))))

      curve(f.lb, add = TRUE, col = sroc.cols[i], lty = sroc.ci.lty, lwd = sroc.ci.lwd)

      curve(f.ub, add = TRUE, col = sroc.cols[i], lty = sroc.ci.lty, lwd = sroc.ci.lwd)
    }

      }


}

