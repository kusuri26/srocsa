###########################
##
## Support function
##
###########################  


sroc.pb.t.fix <- function(
  y1, y2, v1, v2, 
  p,
  c11 = 0.5,
  ci.level = 0.95,
  parset = list()
){

  n  <- length(y1)

    ## fix c11


  c22 <- 1-c11
  c1  <- sqrt(c11)
  c2  <- sqrt(c22)

  ## extract initial values ----
  start5 <- sroc.init(y1, y2, v1, v2, v12=NULL)
  start6 <- c(start5, beta=0.1)

  parset <- modifyList(default.parset.sroc(), parset)


  fn <- function(par) llk.o.srocpb.t(
    par = c(par[1:6], c11),
    n=n, y1=y1, y2=y2, v1=v1, v2=v2, 
    p = p,
    alpha.bound = parset$alpha.bound
    )

  opt <- try(

    nlminb(
    start6,
    fn,
    lower = c(-parset$mu.bound, -parset$mu.bound, parset$eps, parset$eps, -1, parset$beta.interval[1]),
    upper = c( parset$mu.bound,  parset$mu.bound, parset$tau.bound, parset$tau.bound,  1, parset$beta.interval[2])
  )
    , silent = TRUE)


  if(!inherits(opt,"try-error")) {

    ##  output: all pars ----

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

    b   <- opt$par[6]

    t.ldor <- c11*t11 + c22*t22 + 2*c1*c2*t12
    u.ldor <- c1*u1 + c2*u2

    se.ldor2 <- c11*v1+c22*v2
    se.ldor  <- sqrt(se.ldor2)

    sq     <- sqrt(1 + b^2 * (1 + t.ldor/se.ldor2))


    ## hessian ----

    hes <- numDeriv::hessian(fn, opt$par)
    rownames(hes) <- colnames(hes) <- c("mu1", "mu2", "tau1", "tau2", "rho", "beta")

    ## sauc and CI ----
    # hes <- opt$num.hessian

    if(p==1) inv.I.fun.m <- solve(hes[1:5,1:5]) else inv.I.fun.m <- solve(hes)

    opt$var.ml <- inv.I.fun.m

    var.matrix <-  inv.I.fun.m[1:5, 1:5]

    opt$sauc.ci <- sauc(mu1 = u1, mu2 = u2, tau1 = t1, tau2 = t2, rho = r, var.matrix = var.matrix, ci.level = ci.level)
    sauc <- opt$sauc.ci[1]

    ## beta ----
    if(p==1) beta <- NA else beta <- b
    opt$beta <- beta

      # b.se <- suppressWarnings(sqrt(inv.I.fun.m[6,6]))
      # b.lb <- b + qnorm((1-ci.level)/2, lower.tail = TRUE)*b.se
      # b.ub <- b + qnorm((1-ci.level)/2, lower.tail = FALSE)*b.se
      # names(opt$beta.ci) <- c("beta", "beta.lb", "beta.ub")

    ## alpha -----
    if(p==1) a.opt <- NA else {
      a.p <- function(a) { sum(1/ pnorm( (a + b * u.ldor/se.ldor) / sq ), na.rm = TRUE) - n/p }
      a.opt.try <- suppressWarnings(try(uniroot(a.p, interval = c(-parset$alpha.bound, parset$alpha.bound), extendInt =  "yes"), silent = TRUE))
      a.opt <- a.opt.try$root
    }
    opt$alpha <- c(alpha = a.opt)

    ## mu1 CI -----
    u1.se <- suppressWarnings(sqrt(inv.I.fun.m[1,1]))
    u1.lb <- u1 + qnorm((1-ci.level)/2, lower.tail = TRUE)*u1.se
    u1.ub <- u1 + qnorm((1-ci.level)/2, lower.tail = FALSE)*u1.se

    opt$mu1.ci <- c(u1, u1.lb, u1.ub, se, plogis(u1.lb), plogis(u1.ub))
    names(opt$mu1.ci) <- c("mu1", "mu1.lb", "mu1.ub", "se", "se.lb", "se.ub")

    ## mu2 CI ----
    u2.se <- suppressWarnings(sqrt(inv.I.fun.m[2,2]))
    u2.lb <- u2 + qnorm((1-ci.level)/2, lower.tail = TRUE)*u2.se
    u2.ub <- u2 + qnorm((1-ci.level)/2, lower.tail = FALSE)*u2.se

    opt$mu2.ci <- c(u2, u2.lb, u2.ub, sp, plogis(u2.lb), plogis(u2.ub))
    names(opt$mu2.ci) <- c("mu2", "mu2.lb", "mu2.ub", "sp", "sp.lb", "sp.ub")

    ## PAR AND ALL PAR -----
    opt$par.all <- c(u1, u2, t11, t22, t12, c11, c22,  beta, a.opt, sauc, se, sp)
    names(opt$par.all) <- c("mu1", "mu2", "tau1^2", "tau2^2", "tau12", "c1^2", "c2^2", "beta", "alpha", "sauc", "sens", "spec")

    ## LOGIT-DATA  -----
    # opt$l.data <- data
    # opt$call <- this.call

    
} else opt <- NULL

  return(opt)

}



### function of fixing 
sroc.pb.t.unfix <- function(
  y1, y2, v1, v2,
  p,
  ci.level = 0.95,
  parset = list()
){

  

  n  <- length(y1)


    ## unfix c11
    start5 <- sroc.init(y1, y2, v1, v2, v12=NULL)
    start7 <- c(start5, beta=0.1, c11=0.1)

    parset <- modifyList(default.parset.sroc(), parset)

    fn <- function(par) llk.o.srocpb.t(
      par,
      n=n, y1=y1, y2=y2, v1=v1, v2=v2,
      p = p,
      alpha.bound = parset$alpha.bound)

    opt <- try(
      nlminb(start7,
               fn,
               lower = c(-parset$mu.bound, -parset$mu.bound, parset$eps, parset$eps,-1, parset$beta.interval[1], 0),
               upper = c( parset$mu.bound,  parset$mu.bound, parset$tau.bound, parset$tau.bound, 1, parset$beta.interval[2], 1)
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

    b   <- opt$par[6]

    c11 <- opt$par[7]
    c1  <- sqrt(c11)
    c22 <- 1-c11
    c2  <- sqrt(c22)


    u.ldor   <- c1*u1 + c2*u2
    t.ldor   <- c11*t11 + c22*t22 + 2*c1*c2*t12

    se.ldor2 <- c11*v1+c22*v2
    se.ldor  <- sqrt(se.ldor2)

    sq     <- sqrt(1 + b^2 * (1 + t.ldor/se.ldor2))

    ##  hessian ----
    hes <- numDeriv::hessian(fn, opt$par)
    rownames(hes) <- colnames(hes) <- c("mu1", "mu2", "tau1", "tau2", "rho", "beta", "c11")

    ## sauc and CI ----
    if(p==1) inv.I.fun.m <- solve(hes[1:5,1:5]) else inv.I.fun.m <- solve(hes)

    opt$var.ml <- inv.I.fun.m

    var.matrix <-  inv.I.fun.m[1:5, 1:5]

    opt$sauc.ci <- sauc(mu1 = u1, mu2 = u2, tau1 = t1, tau2 = t2, rho = r, var.matrix = var.matrix, ci.level = ci.level)
    sauc <- opt$sauc.ci[1]

    # beta -----
    if(p==1) beta <- NA else beta <- b
    opt$beta <- beta

      # b.se <- suppressWarnings(sqrt(solve(hes)[6,6]))
      # b.lb <- b + qnorm((1-ci.level)/2, lower.tail = TRUE)*b.se
      # b.ub <- b + qnorm((1-ci.level)/2, lower.tail = FALSE)*b.se

      # opt$beta.ci <- c(b, b.lb, b.ub)
    # names(opt$beta.ci) <- c("beta", "beta.lb", "beta.ub")

    ## alpha -----
    if(p==1) a.opt <- NA else {
      a.p <- function(a) { sum(1/ pnorm( (a + b * u.ldor/se.ldor) / sq ), na.rm = TRUE) - n/p }
      a.opt.try <- suppressWarnings(try(uniroot(a.p, interval = c(-parset$alpha.bound, parset$alpha.bound), extendInt = "yes"), silent = TRUE))
      a.opt <- a.opt.try$root
    }

    opt$alpha <- c(alpha = a.opt)

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
    opt$par.all <- c(u1, u2, t11, t22, t12, c11, c22, beta, a.opt, sauc, se, sp)
    names(opt$par.all) <- c("mu1", "mu2", "tau1^2", "tau2^2", "tau12", "c1^2", "c2^2", "beta", "alpha", "sauc", "sens", "spec")
  
  } else opt <- NULL

  
  return(opt)

}



##
## LIKELIHOOD OF THE OBSERVED (NEGATIVE)
## NOT OUTPUT
##

llk.o.srocpb.t <- function(
  par,
  n, y1, y2, v1, v2,
  p,
  alpha.bound
){


  u1 <- par[1]
  u2 <- par[2]
  t1 <- par[3]
  t2 <- par[4]
  r  <- par[5]
  b  <- par[6]
  c11 <- par[7]  ## CAN BE EITHER PAR OR GIVEN VALUE

  t11 <- t1^2
  t22 <- t2^2
  t12 <- t1*t2*r

  c1 <- sqrt(c11)
  c22 <- 1-c11
  c2  <- sqrt(c22)

  ldor     <- c1*y1 + c2*y2
  se.ldor2 <- c11*v1+c22*v2
  se.ldor  <- sqrt(se.ldor2)

  u.ldor   <- c1*u1 + c2*u2
  t.ldor   <- c11*t11 + c22*t22 + 2*c1*c2*t12

  t        <- ldor/se.ldor

  ##
  ## FUNCTOIN b(Sigma) ----
  ##

  f.b <- function(a){

    sq <- suppressWarnings(sqrt(1 + b^2 * (1 + t.ldor/se.ldor2)))

    pnorm( (a + b * u.ldor/se.ldor) / sq )

  }


  ##
  ## FIND THE ROOT OF a = a.opt ----
  ##

  a.p <- function(a) {mean(1/f.b(a), na.rm = TRUE) - 1/p}

  a.opt.try <- suppressWarnings(try(uniroot(a.p, interval=c(-alpha.bound, alpha.bound), extendInt="yes"), silent = TRUE)) 

  a.opt <- a.opt.try$root


  ##
  ##  LOGLIKELIHOOD-1 OF y|Sigma ----
  ##

  det.vec <- (v1+t11)*(v2+t22)-t12^2

  log.det.vec <- suppressWarnings(log(det.vec))

  f.l1  <- ((y1-u1)^2*(v2+t22) - 2*(y2-u2)*(y1-u1)*t12 + (y2-u2)^2*(v1+t11)) / det.vec + log.det.vec

  s.l1  <- -0.5*sum(f.l1, na.rm = TRUE)


  ##
  ##  LOGLIKELIHOOD-2 OF a(a.opt) ----
  ##


  f.l2 <- pnorm(a.opt + b * t)

  s.l2 <- sum( log(f.l2), na.rm = TRUE )


  ##
  ##  LOGLIKELIHOOD-3 OF b(a.opt) ----
  ##

  f.l3 <- f.b(a.opt)

  s.l3 <- sum( log(f.l3), na.rm = TRUE )


  ##
  ##  FINAL LOGLIKELIHOOD ----
  ##

  return(-(s.l1 + s.l2 - s.l3)) ## NEGATIVE

}
