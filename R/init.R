#### input validations

resolve_inputs <- function(
  exprs,
  data = NULL,
  env = parent.frame(),
  numeric = TRUE,
  same.length = TRUE
) {
  if (is.null(names(exprs))) {
    stop("`exprs` must be a named list.")
  }

  eval.env <- if (is.null(data)) env else data

  values <- lapply(names(exprs), function(nm) {
    value <- tryCatch(
      eval(exprs[[nm]], envir = eval.env, enclos = env),
      error = function(e) {
        stop("Could not evaluate `", nm, "`: ", conditionMessage(e))
      }
    )

    if (is.null(value)) {
      stop("`", nm, "` must be supplied.")
    }

    value
  })
  names(values) <- names(exprs)

  if (numeric) {
    not.numeric <- names(values)[!vapply(values, is.numeric, logical(1))]

    if (length(not.numeric) > 0) {
      stop(
        "The following inputs must be numeric: ",
        paste(not.numeric, collapse = ", ")
      )
    }
  }

  if (same.length && length(values) > 1) {
    lengths <- lengths(values)

    if (length(unique(lengths)) != 1) {
      stop(
        "These inputs must have the same length: ",
        paste(names(values), collapse = ", ")
      )
    }
  }

  values
}


## inital values

sroc.init <- function(y1, y2, v1, v2, v12 = NULL){

  n <- length(y1)
  if(is.null(v12)) v12 <- rep(0, n)

    fit.m <- mixmeta::mixmeta(cbind(y1,y2),S=cbind(v1, v12, v2), method="ml")
    if(!inherits(fit.m, "try-error")) {
    if(fit.m$converged){
      p1 <- sqrt(fit.m$Psi[1])
      p2 <- sqrt(fit.m$Psi[4])
      p.r<- fit.m$Psi[3]/(p1*p2)

      ## use the converged estimates from mixmeta as initial values
      start5 <- c(fit.m$coefficients, p1, p2, p.r)

      }} else start5 <- c(rep(0.1, 4), -0.1)

      names(start5) <- c("mu1", "mu2", "tau1", "tau2", "rho")

      return(start5)

 
}



tsroc.pb.init <- function(
  y1, y2, y3,
  v1, v2, v3, 
  v12, v13, v23){

  y <- as.matrix( cbind(y1, y2, y3) )
  S <- as.matrix( cbind(v1, v12, v13, v2, v23, v3) )

  fit.m <- try(mixmeta::mixmeta(y, S, method = "ml"),silent = FALSE)

  if((!inherits(fit.m, "try-error")) && fit.m$converged) {

      u0    <- fit.m$coefficients
      tau   <- c(fit.m$Psi)
      t1230 <- sqrt(tau[c(1, 5, 9)])
      r10   <- prod(tau[2]/prod(sqrt(tau[c(1, 5)])))
      r30   <- prod(tau[3]/prod(sqrt(tau[c(1, 9)])))
      r20   <- prod(tau[6]/prod(sqrt(tau[c(5, 9)])))

      start9 <- c(u0, t1230, r10, r20, r30)
      } else start9 <- c(rep(0.1,6), rep(-0.1,3))

      names(start9) <- c("mu1", "mu2", "mu3", "tau1", "tau2", "tau3","rho1", "rho2", "rho3")

      return(start9)

 
}