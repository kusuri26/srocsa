#'@title
#'Use formula to fit CIFs for time-to-event data with intercurrent events
#'
#'
#'@description
#'This function
#'
#'
#'@param  additional parameters setting:
#'
#'
#'@return A list including the fitted object and input variables.
#'
#'@examples
#'
#'@details
#'
#'@references
#'
#'
#' @seealso \code{\link[tteICE]{surv.boot}}, \code{\link[tteICE]{scr.tteICE}}
#'
#' @importFrom stats model.frame model.response terms as.formula
#'
#'
#'
#'

default.parset.sroc <- function() {
  list(
    inital.values = c(rep(0.1,4),-0.1),  ## initial values for mu1, mu2, tau1, tau2, rho
    mu.bound = 5,
    tau.bound = 5,
    beta.interval = c(0, 5),
    alpha.bound = 5,
    eps = 1e-5
  )
}


default.parset.tsroc.pb <- function() {
  list(
    inital.values = c(rep(0.1, 6), rep(-0.1,3)),  ## initial values for mu1, mu2, tau1, tau2, rho
    mu.bound = 5,
    tau.bound = 5,
    beta.interval = c(0, 5),
    alpha.bound = 5,
    eps = 1e-5
  )
}
