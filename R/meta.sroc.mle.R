###########################
##
## Support function
##
###########################  

meta.sroc.mle <- function(
  # data input
  y1, y2, v1, v2, v12 = NULL, 
  ci.level = 0.95,
  parset = list()
){

  n  = length(y1)
  if (is.null(v12))  v12 <- rep(0, n)

  ## update parameter setting
  parset <- modifyList(default.parset.sroc(), parset)
      

 # ## estimation ----

  if(is.null(parset$inital.values) || length(na.omit(parset$inital.values))!=5) parset$inital.values <- sroc.init(y1, y2, v1, v2, v12)

  fn <- function(par) if(is.null(v12)) llk.p.sroc(par, y1, y2, v1, v2) else llk.p.tsroc(par, y1, y2, v1, v2, v12)
    
    opt <- try(

      nlminb(parset$inital.values,
             fn,
             lower = c(-parset$mu.bound, -parset$mu.bound, parset$eps, parset$eps,-1),
             upper = c( parset$mu.bound,  parset$mu.bound, parset$tau.bound, parset$tau.bound, 1)
    )
      ,silent = TRUE)

    if(!inherits(opt,"try-error")) {

    u1  <- opt$par[1]
    se  <- plogis(u1)
    u2  <- opt$par[2]
    sp  <- plogis(u2)

    t1  <- opt$par[3]
    t11 <- t1^2
    t2  <- opt$par[4]
    t22 <- t2^2
    r   <- opt$par[5]
    t12 <- t1*t2*r


    ##  hessian ----
    hes <- numDeriv::hessian(fn, opt$par)
    rownames(hes) <- colnames(hes) <- c("mu1", "mu2", "tau1", "tau2", "rho")

    ## sauc and CI ----
    inv.I.fun.m <- solve(hes)
    opt$var.ml <- inv.I.fun.m
    var.matrix <-  inv.I.fun.m

    opt$sauc.ci <- sauc(mu1 = u1, mu2 = u2, tau1 = t1, tau2 = t2, rho = r, var.matrix = var.matrix, ci.level = ci.level)
    sauc <- opt$sauc.ci[1]


    ## mu1 CI -----
    u1.se <- suppressWarnings(sqrt(inv.I.fun.m[1,1]))
    u1.lb <- u1 + qnorm((1-ci.level)/2, lower.tail = TRUE)*u1.se
    u1.ub <- u1 + qnorm((1-ci.level)/2, lower.tail = FALSE)*u1.se

    opt$mu1.ci <- c(u1, u1.lb, u1.ub, se, plogis(u1.lb), plogis(u1.ub))
    names(opt$mu1.ci) <- c("mu1", "mu1.lb", "mu1.ub", "sens", "se.lb", "se.ub")

    ## mu2 CI ----
    u2.se <- suppressWarnings(sqrt(inv.I.fun.m[2,2]))
    u2.lb <- u2 + qnorm((1-ci.level)/2, lower.tail = TRUE)*u2.se
    u2.ub <- u2 + qnorm((1-ci.level)/2, lower.tail = FALSE)*u2.se

    opt$mu2.ci <- c(u2, u2.lb, u2.ub, sp, plogis(u2.lb), plogis(u2.ub))
    names(opt$mu2.ci) <- c("mu2", "mu2.lb", "mu2.ub", "spec", "sp.lb", "sp.ub")

    ## PAR AND ALL PAR -----
    opt$par.all <- c(u1, u2, t11, t22, t12, sauc, se, sp)
    names(opt$par.all) <- c("mu1", "mu2", "tau1^2", "tau2^2", "tau12",  "sauc", "sens", "spec")
  }


  return(opt)

}


##
## LIKELIHOOD OF THE OBSERVED (NEGATIVE)
## NOT OUTPUT
##

llk.p.sroc <- function(
  par, y1, y2, v1, v2
){

  u1 <- par[1]
  u2 <- par[2]
  t1 <- par[3]
  t2 <- par[4]
  r  <- par[5]
  

  t11 <- t1^2
  t22 <- t2^2
  t12 <- t1*t2*r


  ##
  ##  LOGLIKELIHOOD-1 OF y|Sigma ----
  ##

  det.vec <- (v1+t11)*(v2+t22)-t12^2

  log.det.vec <- suppressWarnings(log(det.vec))

  f.l1  <- ((y1-u1)^2*(v2+t22) - 2*(y2-u2)*(y1-u1)*t12 + (y2-u2)^2*(v1+t11)) / det.vec + log.det.vec

  s.l1  <- -0.5*sum(f.l1, na.rm = TRUE)



  ##
  ##  FINAL LOGLIKELIHOOD ----
  ##

  return(-(s.l1)) ## NEGATIVE

}


llk.p.tsroc <- function(par,y1, y2, v1, v2, v12){

  u1 <- par[1]
  u2 <- par[2]
  
  t1 <- par[3]
  t2 <- par[4]
  
  r  <- par[5]

  t11 <- t1^2
  t22 <- t2^2
  t12 <- t1*t2*r
  
  v11 <- v1  + t11
  v22 <- v2  + t22
  v12 <- v12 + t12


  ##
  ##  LOGLIKELIHOOD-1 OF y|Sigma ----
  ##

  det.vec <- (v11)*(v22)-v12^2

  log.det.vec <- suppressWarnings(log(det.vec))

  f.l1  <- ((y1-u1)^2*(v22) - 2*(y2-u2)*(y1-u1)*v12 + (y2-u2)^2*(v11)) / det.vec + log.det.vec

  s.l1  <- -0.5*sum(f.l1, na.rm = TRUE)

  -s.l1
}