###########################
##
## Support function
## for print.srocsa, summary.srocsa, plot.srocsa
##
###########################


## type of the srocsa object ----
## "srocsa.t"   : list over p (+ $call)
## "tsroc"      : list over time points (no selection, p = 1)
## "tsrocsa.lrt": list over time points, each a list over p

.srocsa.type <- function(x){

  type <- attr(x, "type")
  if (!is.null(type)) return(type)

  ## objects created before the attribute "type" was added
  x <- unclass(x)
  if (!is.null(x[["call"]])) return("srocsa.t")
  if (all(vapply(x, function(z) is.list(z) && is.null(z$par), logical(1)))) "tsrocsa.lrt" else "tsroc"

}


## split the object into blocks (one block per time point) ----

.srocsa.blocks <- function(x){

  type    <- .srocsa.type(x)
  n.study <- attr(x, "n.study")
  x       <- unclass(x)

  n.at <- function(i) if (is.null(n.study)) NA else n.study[[i]]

  if (type == "srocsa.t") {

    fits <- x[names(x) != "call"]
    p    <- attr(x, "sa.p")
    if (is.null(p)) p <- rep(NA, length(fits))

    return(list(.srocsa.block(fits, p = p, n = n.at(1), time = NULL)))
  }

  lapply(seq_along(x), function(i) {

    fits <- if (type == "tsroc") list(x[[i]]) else x[[i]]
    p    <- if (type == "tsroc") 1 else as.numeric(names(fits))

    .srocsa.block(fits, p = p, n = n.at(i), time = names(x)[i])
  })

}


.srocsa.block <- function(fits, p, n, time){

  est <- do.call(rbind, lapply(fits, .srocsa.fit.est))

  list(time = time, n = n, p = p, est = est, fits = fits)

}


## estimates and SE from one fitted model ----
## par.all and var.ml are used by name, because the order of par differs between
## srocsa.t (mu1, mu2, tau1, tau2, rho, beta, c11),
## meta.sroc.mle (mu1, mu2, tau1, tau2, rho) and
## tsrocsa.lrt.single (mu1-3, tau1-3, rho1-3, beta)

.srocsa.fit.est <- function(fit){

  nm  <- c("conv", "sauc", "sauc.se", "sauc.lb", "sauc.ub",
           "mu1", "mu1.se", "mu2", "mu2.se", "tau1", "tau1.se", "tau2", "tau2.se",
           "rho", "rho.se", "beta", "beta.se", "alpha", "alpha.se")
  out <- setNames(rep(NA_real_, length(nm)), nm)
  out["conv"] <- 0

  ## NULL or try-error: estimation failed
  if (!is.list(fit) || is.null(fit$par.all)) return(out)

  pa <- fit$par.all
  V  <- fit$var.ml

  se <- function(v) if (!is.null(V) && v %in% rownames(V)) suppressWarnings(sqrt(V[v, v])) else NA

  ## 1 indicates converged (nlminb: convergence = 0)
  out["conv"] <- as.numeric(!is.null(fit$convergence) && fit$convergence == 0)

  ## tau: standard deviation
  out["mu1"]  <- pa[["mu1"]]
  out["mu2"]  <- pa[["mu2"]]
  out["tau1"] <- sqrt(pa[["tau1^2"]])
  out["tau2"] <- sqrt(pa[["tau2^2"]])
  out["rho"]  <- pa[["tau12"]] / (out[["tau1"]] * out[["tau2"]])

  rho.nm <- if ("rho" %in% rownames(V)) "rho" else "rho1"
  v5     <- c("mu1", "mu2", "tau1", "tau2", rho.nm)

  out[c("mu1.se", "mu2.se", "tau1.se", "tau2.se", "rho.se")] <- vapply(v5, se, numeric(1))

  if ("beta" %in% names(pa)) {
    out["beta"]    <- pa[["beta"]]
    out["beta.se"] <- se("beta")
  }

  if ("alpha" %in% names(pa)) out["alpha"] <- pa[["alpha"]]

  ## SE of alpha: NOT SUPPORTED YET
  ## alpha is not a free parameter; it is solved from mean(1/P_i) = 1/p given the other
  ## parameters, so it is not included in var.ml. A delta-method SE (implicit function)
  ## needs v1, v2 of each study, which are not kept in the fitted object.
  # out["alpha.se"] <- ...

  ## SAUC and CI ----
  if (!is.null(fit$sauc.ci)) {
    out[c("sauc", "sauc.lb", "sauc.ub")] <- fit$sauc.ci[1:3]
  } else if (!is.null(fit$sauc)) {
    out["sauc"] <- fit$sauc
  }

  ## SE of SAUC by the delta method (same variance as used for sauc.ci)
  if (!is.null(V) && all(v5 %in% rownames(V))) {
    did <- try(.DID.sroc(out[["mu1"]], out[["mu2"]], out[["tau1"]], out[["tau2"]], out[["rho"]], V[v5, v5]), silent = TRUE)
    if (!inherits(did, "try-error")) out["sauc.se"] <- suppressWarnings(sqrt(as.numeric(did)))
  }

  out

}


## print one block ----

.srocsa.print.block <- function(b, digits = 3){

  e   <- b$est
  ok  <- !is.na(e[, "mu1"])
  lab <- paste0("p = ", b$p)

  f <- function(v) ifelse(is.na(v), "NA", formatC(v, digits = digits, format = "f"))

  ## "-": not applicable (p = 1 or the model without selection)
  cell <- function(est, se) ifelse(is.na(est), ifelse(ok, "-", "NA"), paste0(f(est), " (", se, ")"))

  ## model optimization ----
  cat("==================================================\n")
  cat("Model optimization\n\n")
  cat("Number of studies :", b$n, "\n")
  if (!is.null(b$time)) cat("time              :", b$time, "\n")
  cat("Optimizer         :", e[, "conv"], "\n\n")
  cat("* 1 indicates converged\n")
  cat("==================================================\n\n")

  ## SAUC ----
  sauc.tab <- rbind(
    c("Marginal p", "SAUC", "SE", "CI.low", "CI.up"),
    cbind(lab, f(e[, "sauc"]), f(e[, "sauc.se"]), f(e[, "sauc.lb"]), f(e[, "sauc.ub"]))
  )
  cat("SAUC Estimates\n")
  .srocsa.cat.table(sauc.tab)
  cat("\n==================================================\n\n")

  ## parameters ----
  par.tab <- rbind(
    c("Marginal p", "mu_se", "mu_sp", "tau_se", "tau_sp", "rho", "beta", "alpha"),
    cbind(lab,
          cell(e[, "mu1"],   f(e[, "mu1.se"])),
          cell(e[, "mu2"],   f(e[, "mu2.se"])),
          cell(e[, "tau1"],  f(e[, "tau1.se"])),
          cell(e[, "tau2"],  f(e[, "tau2.se"])),
          cell(e[, "rho"],   f(e[, "rho.se"])),
          cell(e[, "beta"],  f(e[, "beta.se"])),
          ## SE of alpha: not supported yet (see .srocsa.fit.est)
          cell(e[, "alpha"], "-"))
  )
  cat("Parameter Estimates (SE)\n")
  .srocsa.cat.table(par.tab)
  cat("\n")

}


## print a character matrix (first row = header) ----

.srocsa.cat.table <- function(tab){

  w    <- apply(tab, 2, function(col) max(nchar(col)))
  rows <- apply(tab, 1, function(r) {
    paste0(sprintf("%-*s", w[1], r[1]), "   ", paste(sprintf("%*s", w[-1], r[-1]), collapse = "   "))
  })
  rule <- strrep("-", max(nchar(rows)))

  cat(paste0(c(rule, rows[1], rule, rows[-1]), "\n"), sep = "")

}
