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
#' @param correct.value Imputation value for ``continuity correction''. The default is 0.5.
#' @param correct.type Two types of ``continuity correction'' (\code{c("single", "all")}). The default is \code{"single"}.
#' \code{"single"}: input value for single study.
#' \code{"all"}: input value for all the cells.
#'
#' @param p Specify the marginal (overall) selection (publication) probability; \eqn{Pr(\text{select}) = p}
#' @param c1.square Pre-specified \eqn{c_1^2}{c1-square}.
#' \eqn{c_2^2 = 1 - c_1^2}{c2-square = 1-c1-square}.
#' @param ci.level The significant value for confidence interval.
#' Area under the Reitsma's SROC curve(\code{"sroc"}) or under the Rutter's HSROC curve (\code{"hsroc"}). See details below.
#'
#' @param estim.set A list of other settings for estimation. Bad initial values and estimation intervals will cause non-convergence results.
#' \itemize{
#'     \item \code{reitsma.par0}: Initial values used for estimating the parameters in the bivariate random effects model (Reitsma's model). It should be changed by a vector of \code{c(mu1, mu2, tau1, tau2, rho)}.
#'     \item \code{beta.interval}: The constraint interval for \eqn{\beta}{beta}. The estimation of \eqn{\beta}{b} will be searched within the interval. Take either the positive interval (>0) or the negative interval (<0).
#'     \item \code{alpha.interval}: The constraint interval for \eqn{\alpha}{alpha}. The root of \eqn{\alpha}{a} will be searched within the interval.
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
#' sa.fit1.fc <- srocpb.p.fc(data = IVD, TP = TP, FN = FN, TN = TN, FP = FP, p = 0.7)
#' sa.fit1.fc
#'
#' sa.fit2.fc <- srocpb.p.fc(ldata = IVD_2, y1 = y1, y2 = y2, v1 = v1, v2 = v2, p = 0.7)
#' sa.fit2.fc
#'
#' ## Specifying data and variables are required.
#' ## The followings are wrong and cause error
#' \dontrun{
#' sa.fit1.fc.wrong <- srocpb.p.fc(data = IVD, TP, FN, TN, FP, p = 0.7)
#'
#' sa.fit2.fc.wrong <- srocpb.p.fc(ldata = IVD_2, y1, y2, v1, v2, p = 0.7)
#' }
#'
#' sa.fit2 <- srocpb.p.fc(
#'   data = IVD, TP = TP, FN = FN, TN = TN, FP = FP, p = 0.7,
#'   correct.type = "all"
#' )
#' sa.fit2
#'
#' sa.fit3 <- srocpb.p.fc(IVD, TP = TP, FN = FN, TN = TN, FP = FP, p = 0.7)
#' sa.fit3
#'
#' @seealso \code{\link[stats]{uniroot}}.
#'

sroc <- function(
  TP, FN, FP, TN,
  data = NULL,
  cc.value = 0.5,
  cc.type = c("single", "all"),
  ci.level = 0.95,
  parset = list()
) {
  cc.type <- match.arg(cc.type)
  inputs <- resolve_inputs(
    exprs = list(
      TP  = substitute(TP),
      FN  = substitute(FN),
      FP  = substitute(FP),
      TN  = substitute(TN)
    ),
    data = data,
    env = parent.frame(),
    numeric = FALSE,
    same.length = TRUE
  )


  TP <- inputs$TP
  FN <- inputs$FN
  FP <- inputs$FP
  TN <- inputs$TN


  ## update parameter setting
  parset <- modifyList(default.parset.sroc(), parset)

  data.conv <- convert.dt.sroc(TP, FN, FP, TN, data = NULL, cc.value = cc.value, cc.type = cc.type)

  meta.sroc.mle(data.conv$y1, data.conv$y2, data.conv$v1, data.conv$v2, ci.level = ci.level, parset = parset)
}
