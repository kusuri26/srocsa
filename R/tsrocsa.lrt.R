#' @title Sensitivity analysis of diagnostic meta-analysis with prespecified contrast vector
#'
#' @description Sensitivity analysis with pre-specified (fixed) contrast vector \eqn{\boldsymbol{c} = (c_1, c_2)}
#'
#' @param data Data with \{TP, FN, TN, FP\} data.
#' @param TP Specify the column of true positives (TP).
#' @param FN Specify the column of false negatives (FN).
#' @param TN Specify the column of true negatives (TN).
#' @param FP Specify the column of false positives (FP).
#'
#' @param ldata Data with logit-transformed sensitivity and specificity data and their corresponding variances.
#' @param y1 Specify the column of logit-transformed sensitivity.
#' @param y2 Specify the column of logit-transformed specificity.
#' @param v1 Specify the column of the variance of logit-transformed sensitivity.
#' @param v2 Specify the column of the variance of logit-transformed specificity.
#'
#' @param cc.value Imputation value for ``continuity correction''. The default is 0.5.
#' @param cc.type Two types of ``continuity correction'' (\code{c("single", "all")}). The default is \code{"single"}.
#' \code{"single"}: input value for single study.
#' \code{"all"}: input value for all the cells.
#'
#' @param p Specify the marginal (overall) selection (publication) probability; \eqn{Pr(\text{select}) = p}
#' @param c1.square Pre-specified \eqn{c_1^2}{c1-square}.
#' \eqn{c_2^2 = 1 - c_1^2}{c2-square = 1-c1-square}.
#' @param ci.level The significant value for confidence interval.
#' Area under the Reitsma's SROC curve(\code{"sroc"}) or under the Rutter's HSROC curve (\code{"hsroc"}). See details below.
#'
#' @param parset A list of other settings for estimation. Bad initial values and estimation intervals will cause non-convergence results.
#' \itemize{
#'     \item \code{reitsma.par0}: Initial values used for estimating the parameters in the bivariate random effects model (Reitsma's model). It should be changed by a vector of \code{c(mu1, mu2, tau1, tau2, rho)}.
#'     \item \code{beta.interval}: The constraint interval for \eqn{\beta}{beta}. The estimation of \eqn{\beta}{b} will be searched within the interval. Take either the positive interval (>0) or the negative interval (<0).
#'     \item \code{alpha.bound}: The constraint interval for \eqn{\alpha}{alpha}. The root of \eqn{\alpha}{a} will be searched within the interval.
#'     \item \code{eps}: A minimum positive value.
#'     \item \code{show.warn.message}: Whether to show the warning messages.
#'   }
#' @param ... See other augments in function \code{\link[stats]{uniroot}}.
#'
#' @return
#' publication bias adjusted estimations,
#' confidence interval,
#' convergence list,
#' logit transformed data
#' called function
#'
#' @details
#' \describe{
#' \item{Continuity correction}{\url{https://en.wikipedia.org/wiki/Continuity_correction}}
#'
#' \item{Reitsma's SROC curve and the Rutter's HSROC curve}
#' }
#'
#' @examples
#'
#' sa.fit1.fc <- srocpb.p.fc(data=IVD, TP = TP, FN = FN, TN = TN, FP = FP, p = 0.7)
#' sa.fit1.fc
#'
#' sa.fit2.fc <- srocpb.p.fc(ldata=IVD_2, y1 = y1, y2 = y2, v1 = v1, v2 = v2, p = 0.7)
#' sa.fit2.fc
#'
#' ## Specifying data and variables are required.
#' ## The followings are wrong and cause error
#' \dontrun{
#' sa.fit1.fc.wrong <- srocpb.p.fc(data=IVD, TP, FN, TN, FP, p = 0.7)
#'
#' sa.fit2.fc.wrong <- srocpb.p.fc(ldata=IVD_2, y1, y2, v1, v2, p = 0.7)
#' }
#'
#' sa.fit2 <- srocpb.p.fc(data=IVD, TP = TP, FN = FN, TN = TN, FP = FP, p = 0.7,
#'                           cc.type = "all")
#' sa.fit2
#'
#' sa.fit3 <- srocpb.p.fc(IVD, TP = TP, FN = FN, TN = TN, FP = FP, p = 0.7)
#' sa.fit3
#'
#' @seealso \code{\link[stats]{uniroot}}.
#'
#' @export

tsrocsa.lrt <- function(
  time.points, 
  study.os, time.os, n1.os, n0.os, s1.os, s0.os, data.os,
  study.lnhr, lnhr, var.lnHR, data.lnhr,
  time.med, n1.med, n0.med, s1.med, s0.med, data.mct,
  senp.p = 1,
  ci.level = 0.95,
  parallel.ncores = 1L,
  parset = list()
){

  #### data.os
  inputs <- resolve_inputs(
      exprs = list(
        study.os  = substitute(study.os),
        time.os     = substitute(time.os),
        n1.os     = substitute(n1.os),
        n0.os     = substitute(n0.os),
        s1.os     = substitute(s1.os),
        s0.os     = substitute(s0.os)
      ),
      data = data.os,
      env = parent.frame(),
      numeric = FALSE,
      same.length = TRUE
    )

    data.os <- data.frame(
      study.os  = inputs$study.os,
      time.os = inputs$time.os,
      n1.os   = inputs$n1.os,
      n0.os   = inputs$n0.os,
      s1.os   = inputs$s1.os,
      s0.os   = inputs$s0.os
    )

#### data.lnhr

  inputs <- resolve_inputs(
      exprs = list(
        study.lnhr  = substitute(study.lnhr),
        lnhr = substitute(lnhr),
        var.lnHR = substitute(var.lnHR)
      ),
      data = data.lnhr,
      env = parent.frame(),
      numeric = FALSE,
      same.length = TRUE
    )

    data.lnhr <- data.frame(
      study.lnhr  = inputs$study.lnhr,
      lnhr = inputs$lnhr,
      var.lnHR = inputs$var.lnHR
    )

#### data.mtc
    inputs <- resolve_inputs(
      exprs = list(
        time.med = substitute(time.med),
        n1.med   = substitute(n1.med),
        n0.med   = substitute(n0.med),
        s1.med   = substitute(s1.med),
        s0.med   = substitute(s0.med)
      ),
      data = data.mct,
      env = parent.frame(),
      numeric = FALSE,
      same.length = TRUE
    )

    data.mct <- data.frame(
      time.med  = inputs$time.med,
      n1.med    = inputs$n1.med,
      n0.med    = inputs$n0.med,
      s1.med    = inputs$s1.med,
      s0.med    = inputs$s0.med
    )


        ## update parameter setting
  parset <- modifyList(default.parset.tsroc.pb(), parset)


f <- function(tt, pp){

  data.log <- convert.dt.tsroc(
    time.point = tt, study.os=study.os, time.os, n1.os, n0.os, s1.os, s0.os, data.os,
    study.lnhr, lnhr, var.lnHR, data.lnhr,
    time.med, n1.med, n0.med, s1.med, s0.med, data.mct)

  tsrocsa.lrt.single(
    y1=data.log$u_sen, y2=data.log$u_spe, y3=data.log$u_lnHR,
    v1=data.log$v_sen, v2=data.log$v_spe, v3=data.log$v_lnHR,
    v12=data.log$v_senspe, v13 = data.log$v_senlnHR, v23=data.log$v_spelnHR,
    p=pp, ci.level=ci.level,
    parset=parset)
  
  }

 ## estimation ----

 res <- NULL

  
  fit_at_time <- function(tt) {
    lapply(senp.p, function(pp) f(tt = tt, pp = pp))
  }

  if (parallel.ncores == 1L) {
    res <- lapply(time.points, fit_at_time)
  } else {
    cl <- parallel::makeCluster(parallel.ncores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    ## PSOCK workers start with empty global environments, so export every
    ## helper reached by the fitting closure.
    parallel::clusterExport(
      cl,
      varlist = c(
        "f", "fit_at_time", "convert.dt.tsroc", "cens.eta",
        "tsrocsa.lrt.single", "tsroc.pb.init", "clk.TNM.ml",
        "meta.sroc.mle", "llk.p.sroc", "llk.p.tsroc", "sroc.init",
        "sauc", ".DID.sroc", "resolve_inputs",
        "default.parset.tsroc.pb", "default.parset.sroc"
      ),
      envir = environment()
    )

    res <- parallel::parLapply(cl, X = time.points, fun = fit_at_time)
  }

  class(res) <- "srocsa"
  names(res) <- as.character(time.points)
  for (i in seq_along(res)) {
    names(res[[i]]) <- as.character(senp.p)
  }

  ## information for print/summary/plot methods
  attr(res, "type")    <- "tsrocsa.lrt"
  attr(res, "n.study") <- vapply(time.points, function(tt) sum(data.os$time.os == tt, na.rm = TRUE), numeric(1))

  return(res)



}





