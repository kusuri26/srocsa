###########################
##
## Support function
##
###########################  

tsrocsa.lrt.single <- function(
  y1, y2, y3,
  v1, v2, v3, v12, v13, v23,
  p, 
  ci.level = 0.95,
  parset = list()
  ){


  ## update parameter setting
  parset <- modifyList(default.parset.tsroc.pb(), parset)
  
  

  if(p==1){
      
      opt.o <- meta.sroc.mle(y1=y1, y2=y2, v1=v1, v2=v2, v12=v12, parset=parset)

    }

    else{


      fn.o <- function(par) clk.TNM.ml(par, y1, y2, y3, v1, v2, v3, v12, v13, v23,
                                        p = p,
                                        a.interval= c(-parset$alpha.bound, parset$alpha.bound))

      if(is.null(parset$inital.values) | length(na.omit(parset$inital.values)) ) {
        start9 <- tsroc.pb.init(y1, y2, y3, v1, v2, v3, v12, v13, v23)
        start10 <- c(start9, 0.1)

      }

      opt.o = nlminb(start10, 
                 fn.o,
                 lower = c(rep(-parset$mu.bound,3),rep(0,3),                rep(-1,3), parset$beta.interval[1]),
                 upper = c(rep( parset$mu.bound,3),rep(parset$tau.bound,3), rep(1,3), parset$beta.interval[2]))

      if(!inherits(opt.o, "try-error")) {

        u1  <- opt.o$par[1]
        se  <- plogis(u1)
        u2  <- opt.o$par[2]
        sp  <- plogis(u2)

        u3  <- opt.o$par[3]
        hr  <- exp(u2)

        t1  <- opt.o$par[4]
        t11 <- t1^2
        t2  <- opt.o$par[5]
        t22 <- t2^2
        t3  <- opt.o$par[6]
        t33 <- t3^2

        r1   <- opt.o$par[7]
        t12  <- t1*t2*r1
        r2   <- opt.o$par[8]
        t23  <- t2*t3*r2
        r3   <- opt.o$par[9]
        t13  <- t1*t3*r3

        b    <- opt.o$par[10]


        ##  hessian ----
        num.hessian = numDeriv::hessian(fn.o, opt.o$par)
        rownames(num.hessian) = colnames(num.hessian) = c(paste0("mu",1:3), paste0("tau",1:3), paste0("rho",1:3), "beta")

        var.ml = solve(num.hessian)
        var.matrix = var.ml[c(1,2,4,5,7), c(1,2,4,5,7)]
        sauc.ci = sauc(mu1 = u1, mu2 = u2, tau1 = t1, tau2 = t2, rho = r1, var.matrix = var.matrix, ci.level = ci.level)
        opt.o$sauc <- sauc <- sauc.ci[1]


        opt.o$beta <- b

        f.b = function(a){
          
          sq = suppressWarnings(sqrt(1 + b^2 * (1 + t33 / v3)))
          pnorm( (a + b * u3/sqrt(v3)) / sq )
          
        }
        
        a.p = function(a) {mean(1/f.b(a), na.rm = TRUE) - 1/p}
        
        a.opt.try = suppressWarnings(try(uniroot(a.p, interval=c(-parset$alpha.bound, parset$alpha.bound), extendInt="yes"), silent = TRUE)) 
        a.opt = a.opt.try$root
        opt.o$alpha <- c(alpha = a.opt)

        ## mu1 CI -----
        u1.se <- suppressWarnings(sqrt(var.ml[1,1]))
        u1.lb <- u1 + qnorm((1-ci.level)/2, lower.tail = TRUE)*u1.se
        u1.ub <- u1 + qnorm((1-ci.level)/2, lower.tail = FALSE)*u1.se
        opt.o$mu1.ci <- c(u1, u1.lb, u1.ub, se, plogis(u1.lb), plogis(u1.ub))
        names(opt.o$mu1.ci) <- c("mu1", "mu1.lb", "mu1.ub", "sens", "se.lb", "se.ub")

        ## mu2 CI ----
        u2.se <- suppressWarnings(sqrt(var.ml[2,2]))
        u2.lb <- u2 + qnorm((1-ci.level)/2, lower.tail = TRUE)*u2.se
        u2.ub <- u2 + qnorm((1-ci.level)/2, lower.tail = FALSE)*u2.se
        opt.o$mu2.ci <- c(u2, u2.lb, u2.ub, sp, plogis(u2.lb), plogis(u2.ub))
        names(opt.o$mu2.ci) <- c("mu2", "mu2.lb", "mu2.ub", "spec", "sp.lb", "sp.ub")

        ## mu3 CI ----
        u3.se <- suppressWarnings(sqrt(var.ml[3,3]))
        u3.lb <- u2 + qnorm((1-ci.level)/2, lower.tail = TRUE)*u3.se
        u3.ub <- u2 + qnorm((1-ci.level)/2, lower.tail = FALSE)*u3.se
        opt.o$mu3.ci <- c(u3, u3.lb, u3.ub, hr, exp(u2.lb), exp(u2.ub))
        names(opt.o$mu2.ci) <- c("mu2", "mu2.lb", "mu2.ub", "hr", "hr.lb", "hr.ub")


        ## PAR AND ALL PAR -----
        opt.o$par.all <- c(u1, u2, u3,
          t11, t22, t33,
          t12, t23, t13,
          b, a.opt, sauc, 
          se, sp, hr)
        names(opt.o$par.all) <- c("mu1", "mu2", "mu3", "tau1^2", "tau2^2", "tau3^2", "tau12", "tau23", "tau13", "beta", "alpha", "sauc", "sens", "spec", "hr")

        ## for print/summary/plot methods (appended at the end of the list)
        opt.o$sauc.ci <- sauc.ci
        opt.o$var.ml  <- var.ml

      } else opt.o <- NULL
    }

    return(opt.o)

  
  }



clk.TNM.ml <- function(
  par, y1, y2, y3,
  v1, v2, v3, v12, v13, v23,
  p, a.interval
){

  u1  <- par[1]      ## logit-se
  u2  <- par[2]      ## logit-sp
  u3  <- par[3]      ## ln-HR
  
  t1  <- par[4]      ## logit-se
  t2  <- par[5]    
  t3  <- par[6]
  
  r1  <- par[7]      ## logit-sp
  r2  <- par[8]
  r3  <- par[9]      ## ln-HR
  
  b   <- par[10] 

  t11 <- t1^2
  t22 <- t2^2
  t33 <- t3^2

  t12 <- t1*t2*r1
  t23 <- t2*t3*r2
  t13 <- t1*t3*r3
  
  v11 <- v1  + t11
  v22 <- v2  + t22
  v33 <- v3  + t33
  
  v12 <- v12 + t12 
  v13 <- v13 + t13
  v23 <- v23 + t23

  t_lnHR   <- y3/sqrt(v3)

  ##
  ## FUNCTOIN b(Sigma) ----
  ##

  f.b <- function(a){
    
    sq <- suppressWarnings(sqrt(1 + b^2 * (1 + t33 / v3)))
    
    pnorm( (a + b * u3/sqrt(v3)) / sq )
    
  }


  ##
  ## FIND THE ROOT OF a = a.opt ----
  ##

  a.p <- function(a) {mean(1/f.b(a), na.rm = TRUE) - 1/p}

  a.opt.try <- suppressWarnings(try(uniroot(a.p, a.interval, extendInt="yes"), silent = TRUE)) 

  a.opt <- a.opt.try$root


  ##
  ##  LOGLIKELIHOOD-1 OF y|Sigma ----
  ##
  
  det1  <- v22 * v33 - v23 * v23
  det2  <- v12 * v33 - v13 * v23
  det3  <- v12 * v23 - v13 * v22
  
  det   <- v11 * det1 - v12 * det2 + v13 * det3
  
  log.dev <- suppressWarnings(log(det)) 
  
  a11 <-  (v22 * v33 - v23 * v23)
  a12 <- -(v12 * v33 - v13 * v23)
  a13 <-  (v12 * v23 - v13 * v22)
  
  a21 <- -(v12 * v33 - v23 * v13)
  a22 <-  (v11 * v33 - v13 * v13)
  a23 <- -(v11 * v23 - v12 * v13)
  
  a31 <-  (v12 * v23 - v13 * v22)
  a32 <- -(v11 * v23 - v12 * v13)
  a33 <-  (v11 * v22 - v12 * v12)
  
  Y1  <- y1 - u1
  Y2  <- y2 - u2
  Y3  <- y3 - u3
  
  y_Adj_y <- Y1 * (Y1 * a11 + Y2 * a21 + Y3 * a31) + 
    Y2 * (Y1 * a12 + Y2 * a22 + Y3 * a32) +
    Y3 * (Y1 * a13 + Y2 * a23 + Y3 * a33)
  
  
  
  s.l1 <- -0.5 * (sum(log.dev + y_Adj_y / det, na.rm = TRUE))
  

  ##
  ##  LOGLIKELIHOOD-2 OF a(a.opt) ----
  ##


  f.l2 <- pnorm(a.opt + b * t_lnHR)

  s.l2 <- sum( log(f.l2), na.rm = TRUE )


  ##
  ##  LOGLIKELIHOOD-3 OF b(a.opt) ----
  ##

  f.l3 <- f.b(a.opt)

  s.l3 <- sum( log(f.l3), na.rm = TRUE )


  ##
  ##  FINAL LOGLIKELIHOOD ----
  ##

  -(s.l1 + s.l2 - s.l3) ## NEGATIVE

  
}
