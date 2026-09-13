#' tteICE: Treatment Effect Estimation for Time-to-Event Data with Intercurrent Events
#'
#' This package aims to analyze treatment effects in clinical trials with time-to-event outcomes is complicated by intercurrent events.
#' This package implements methods for estimating and inferring the cumulative incidence functions for time-to-event (TTE) outcomes with intercurrent events (ICE) under the five strategies outlined in the ICH E9 (R1) addendum, see Deng (2025) <doi:10.1002/sim.70091>.
#' This package can be used for analyzing data from both randomized controlled trials and observational studies.
#' In general, the data involve a primary outcome event and, potentially, an intercurrent event.
#' Two data structures are allowed: competing risks, where only the time to the first event is recorded,
#' and semicompeting risks, where the times to both the primary outcome event and intercurrent event (or censoring) are recorded.
#' For estimation methods, nonparametric estimation (which does not use covariates) and semiparametrically efficient estimation are presented.
#'
#' Main functions:
#' \itemize{
#' 	 \item \code{\link[tteICE]{tteICE}} Using formula to fit cumulative incidence functions (CIFs) for competing/semicompeting risk time-to-event data with intercurrent events.
#'   \item \code{\link[tteICE]{scr.tteICE}} Fit CIFs for semicompeting risk time-to-event data with intercurrent events.
#'   \item \code{\link[tteICE]{surv.tteICE}} Fit CIFs for competing risk time-to-event with intercurrent events.
#'   \item \code{\link[tteICE]{tteICEShiny}} Interactive Shiny app for the 'tteICE' package
#'   }
#'
#' Results output functions:
#' \itemize{
#' 	 \item \code{\link[tteICE]{plot.tteICE}} Plot results from 'tteICE' objects.
#'   \item \code{\link[tteICE]{print.tteICE}} Print a short summary of results from 'tteICE' objects
#'   \item \code{\link[tteICE]{summary.tteICE}} Summarize results from 'tteICE' objects
#'   \item \code{\link[tteICE]{predict.tteICE}} Predict risks for 'tteICE' objects at specific time 
#'   \item \code{\link[tteICE]{coef.tteICE}} Show the coefficient for 'tteICE' objects
#'   \item \code{\link[tteICE]{bshaz.tteICE}} Extract the baseline hazards for 'tteICE' objects
#'   \item \code{\link[tteICE]{zph.tteICE}} Perform a test for the proportional hazards assumption for the Cox models in 'tteICE' objects
#'
#'
#' }
#' Example data:
#' \itemize{
#'   \item \code{\link{bmt}} Data from Section 1.3 of Klein and Moeschberger (1997)
#' }
#'
#' @docType package
#' @name tteICE-package
#'
#' @importFrom lifecycle deprecated
"_PACKAGE"
