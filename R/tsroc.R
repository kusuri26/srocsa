#' @title Time Dependent SROC without PB
#'
#' @description Time Dependent Sensitivity Analysis Old (not recommended)
#'
#' @param HR.data Data with variable names
#' \{study, HR, ci.low, ci.up, n1, n0\} .
#' If not, please Set the name in a later argument or change the variable names.
#' 'study' should be a number that is unified with the two data below
#'
#' @param OS.data Data with variable names
#' \{study, t, n1, n0, s1, s0\} .
#' If not, please Set the name in a later argument or change the variable names.
#' \{n1,n0\} is same with HR.data.Either can be omitted.
#'
#' @param MCT.data Data with variable names
#' \{study, n1.mct, n0.mct, s1.med.mct, s0.med.mct, med.year\}
#' If not, please Set the name in a later argument or change the variable names.
#'
#'
#' @param p Specified probability of selection (or publication); Pr(select) = p
#'
#' @param period Convert the units of t to years. When the unit is month, use 12 as an argument.
#'
#' @param ci.level The significant value for confidence interval of SAUC.
#'
#' @param eta.range eta range
#'
#' @param init.eta initial eta
#'
#' @param sauc.type Two types of SAUC values.
#' Area under the Reitsma's SROC curve(\code{"sroc"}) or under the Rutter's HSROC curve (\code{"hsroc"}).
#'
#'
#' @param study  The three input data sets are unified with the same \code{study} number.
#' Use the name of the data field as an argument
#'
#' @param time time column name
#'
#' @param n1  OS.data n1 name
#' @param n0  OS.data n0 name
#' @param s1  OS.data s1 name
#' @param s0  OS.data s0 name
#'
#' @param n1.mct1  MCT.data n1 name
#' @param n0.mct  MCT.data n0 name
#' @param s1.med.mct  MCT.data s1 name
#' @param s0.med.mct  MCT.data s0 name
#'
#' @param med.year MCT.data time
#'
#' @param ... See other augments in function \code{\link[stats]{uniroot}}.
#'
#' @return
#' mu tau
#'
#' @examples
#'
#' res<-tdsaold(HR,OS,MCT,p=c(1,0.8,0.6),s1.med.mct = s1_mct,s0.med.mct = s0_mct,period = 12,med.year = mct_mo)
#' plot(res)
#'
#' @details
#' \describe{
#' \item{Continuity correction}{\url{https://en.wikipedia.org/wiki/Continuity_correction}}
#' }
#'
#' @seealso \code{\link[stats]{uniroot}},
#' \code{\link{plot.tdsaold}}.
#'
#' 
tsroc <- function(
  time.points, 
  study.os, time.os, n1.os, n0.os, s1.os, s0.os, data.os=NULL,
  study.lnhr, lnhr, var.lnHR, data.lnhr =NULL,
  time.med, n1.med, n0.med, s1.med, s0.med, data.mct=NULL,
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


  f <- function(tt){

  data.log <- convert.dt.tsroc(
    time.point = tt, study.os, time.os, n1.os, n0.os, s1.os, s0.os, data.os,
    study.lnhr, lnhr, var.lnHR, data.lnhr,
    time.med, n1.med, n0.med, s1.med, s0.med, data.mct)

  meta.sroc.mle(y1=data.log$u_sen, y2=data.log$u_spe, v1=data.log$v_sen, v2=data.log$v_spe, v12=data.log$v_senspe, ci.level = ci.level, parset=parset)
  
  }

  if(parallel.ncores == 1L) {

     res <- lapply(time.points, function(tt) f(tt))

  } else {

    cl <- parallel::makeCluster(parallel.ncores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      varlist = c("f", "convert.dt.tsroc", "meta.sroc.mle", "sroc.init"),
      envir = environment()
    )

    ## These functions were defined outside metacont.pb(), usually in .GlobalEnv
    parallel::clusterExport(
      cl,
      varlist = c("resolve_inputs", "default.parset.sroc"),
      envir = .GlobalEnv
    )

    res <- parallel::parLapply(
      cl,
      X = time.points,
      fun = f
    )


  }

# res$call <- this.call
  
  class(res) <- "srocsa"
  names(res) <- time.points

  ## information for print/summary/plot methods
  attr(res, "type")    <- "tsroc"
  attr(res, "n.study") <- vapply(time.points, function(tt) sum(data.os$time.os == tt, na.rm = TRUE), numeric(1))

  return(res)

# setClass(

#   Class = "tdsameta",
  
#   slot=c(data="data.frame",par="matrix",var.matrix="list")

# )

}