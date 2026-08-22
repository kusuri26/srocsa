#' @title Print dtametasa.fc and dtametasa.rc functions results
#'
#' @description Print results from function \code{\link{dtametasa.fc}} or \code{\link{dtametasa.rc}}
#'
#' @param x object from function \code{dtametasa.fc} or \code{dtametasa.rc}
#' @param digits digits of the results
#' @param ... other parameters in function \code{\link{print}}
#'
#' @importFrom stats integrate nlminb plogis pnorm qlogis uniroot qchisq qnorm
#' @importFrom mixmeta mixmeta
#'
#'
#' @seealso
#' \code{\link{dtametasa.fc}},
#' \code{\link{dtametasa.rc}},
#' \code{\link[base]{print}}.
#'
#' @rdname print.dtametasa
#'
#' @export

print.meta.pb <- function(x, digits = 3, ...){


  print(list(par.all = x$par.all), digits = digits, ...)

}
