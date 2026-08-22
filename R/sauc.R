#' @title Summary ROC curves
#'
#' @description caluculate SAUC and Confidence interval
#'
#' @param ... (1) The result from function \code{cover.dt}or\code{\link[dtametasa]{dtametasa.fc}} or \code{\link[dtametasa]{dtametasa.rc}} (single curve);
#' (2) a matrix of which the rows are \code{c(mu1, mu2, tau1, tau2, rho)} (multiple curves);
#' (3) a vector of \code{c(mu1, mu2, tau1, tau2, rho)} (single curve). 
#' (4) as arguments (for example\code{mu1=0.8,mu2=0.7,tau1=0.8,tau2=0.9,rho=-0.3} )argument name is important.
#'
#' @param var.matrix variance matrix
#' 
#' @param sroc.type Plot Reitsma's SROC curve(\code{"sroc"}) or Rutter's HSROC curve (\code{"hsroc"}),
#' 
#' @param ci.level Significant level of the confidence intervals.
#'
#'
#' @return Sauc and CI
#'
#' @seealso
#' \code{\link[dtametasa]{dtametasa.fc}},
#' \code{\link{tdsameta}}.
#'
#'
#' @examples
#' res<-tdsameta(HR,OS,MCT,s1.med.mct = s1_mct,s0.med.mct = s0_mct,period = 12,med.year = mct_mo,prob=c(1,0.8))
#' res@@par[c("sauc.sauc", "sauc.sauc.lb", "sauc.sauc.ub"),]
#' 
#'
#' @export


sauc <- function(
	mu1, mu2, tau1, tau2, rho,
	var.matrix = NULL,
	ci.level){


		sroc <- function(x) plogis(mu1 - (tau1*tau2*rho/(tau2^2)) * (qlogis(x) + mu2))
		sauc.try <- try(integrate(sroc, 0, 1))

		if(!inherits(sauc.try, "try-error")) sauc <- sauc.try$value else sauc <- NA

		if(!is.null(var.matrix)){
		sauc.lb <-  plogis(qlogis(sauc) + qnorm((1-ci.level)/2, lower.tail = TRUE) *
																							suppressWarnings(
																								sqrt(.DID.sroc(mu1, mu2, tau1, tau2, rho, var.matrix))/(sauc*(1-sauc))) )

		sauc.ub <-  plogis(qlogis(sauc) + qnorm((1-ci.level)/2, lower.tail = FALSE)*
																							suppressWarnings(
																								sqrt(.DID.sroc(mu1, mu2, tau1, tau2, rho, var.matrix))/(sauc*(1-sauc))) )


		sauc.ci <- c(sauc, sauc.lb, sauc.ub)
		names(sauc.ci) <- c("sauc", "sauc.lb", "sauc.ub")

		} else {

			sauc.ci <- c(sauc, NA, NA)
			names(sauc.ci) <- c("sauc", "sauc.lb", "sauc.ub")
		}

sauc.ci
}




# ##
# ## DID for HSROC (not change yet)
# ##

# .DID.hsroc <- function(mu1, mu2, tau1, tau2, rho, var.matrix){

# 	Q1 <- function(x) {

# 		g <- plogis(mu1 + mu2*tau1/tau2 + tau1/tau2*qlogis(x))

# 		p.mu1 <- g*(1-g) * 1
# 	}

# 	Q2 <- function(x) {

# 		g <- plogis(mu1+mu2*tau1/tau2 + tau1/tau2*qlogis(x))

# 		p.mu2 <- g*(1-g) * (tau1/tau2)
# 	}

# 	Q3 <- function(x) {

# 		g <- plogis(mu1+mu2*tau1/tau2 + tau1/tau2*qlogis(x))
# 		p.tau1 <- g*(1-g) * (mu2 + qlogis(x))/tau2

# 	}

# 	Q4 <- function(x) {

# 		g <- plogis(mu1+mu2*tau1/tau2 + tau1/tau2*qlogis(x))
# 		p.tau2 <- g*(1-g) * (-1/tau2^2)*(mu2*tau1+tau1*qlogis(x))

# 	}



# 	fd <- c(integrate(Q1, 0, 1)$value,
# 									integrate(Q2, 0, 1)$value,
# 									integrate(Q3, 0, 1)$value,
# 									integrate(Q4, 0, 1)$value,
# 									0
# 	)

# 	(fd %*% var.matrix %*% fd)

# }



##
## DID for SROC
##

.DID.sroc <- function(mu1, mu2, tau1, tau2, rho, var.matrix){

	Q1 <- function(x) {

		g <- plogis(mu1 - (tau1*tau2*rho/(tau2^2)) * (qlogis(x) + mu2))
		g*(1-g)
	}

	Q2 <- function(x) {

		g <- plogis(mu1 - (tau1*tau2*rho/(tau2^2)) * (qlogis(x) + mu2))
		p.mu2 <- (-rho*tau1/tau2)*g*(1-g)
	}

	Q3 <- function(x) {

		g <- plogis(mu1 - (tau1*tau2*rho/(tau2^2)) * (qlogis(x) + mu2))
		p.tau1 <- (-rho/tau2*(qlogis(x)+mu2))*g*(1-g)

	}

	Q4 <- function(x) {

		g <- plogis(mu1 - (tau1*tau2*rho/(tau2^2)) * (qlogis(x) + mu2))
		p.tau2 <- rho*tau1/tau2^2*( qlogis(x)+mu2)*g*(1-g)

	}

	Q5 <- function(x) {

		g <- plogis(mu1 - (tau1*tau2*rho/(tau2^2)) * (qlogis(x) + mu2))
		p.rho  <- (-tau1)/tau2*(qlogis(x) + mu2)*g*(1-g)
	}

	fd <- c(integrate(Q1, 0, 1)$value,
									integrate(Q2, 0, 1)$value,
									integrate(Q3, 0, 1)$value,
									integrate(Q4, 0, 1)$value,
									integrate(Q5, 0, 1)$value
	)

	(fd %*% var.matrix %*% fd)

}



##
## FOR SROC CI
##

.QIQ.sroc <- function(x, mu1, mu2, tau1, tau2, rho, inv.I.fun.m) {

  sapply(1: length(x), function(i){

    Q <- c(1, -rho*tau1/tau2, -rho/tau2*(qlogis(x[i])+mu2), rho*tau1/tau2^2*(qlogis(x[i])+mu2), -tau1/tau2*(qlogis(x[i])+mu2))

    (Q %*% inv.I.fun.m %*% Q)

  })

}



