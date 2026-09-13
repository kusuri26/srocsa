###########################
##
## Support function
##
###########################  
convert.dt.sroc <- function(
  TP, FN, FP, TN, data,
  cc.value = 0.5,
  cc.type  = c("single", "all")
  ){

  cc.type <- match.arg(cc.type)

    if (cc.value <=0 | cc.value>1) {
    stop("The value of correct.value must be positive and within (0, 1]",  call. = FALSE)
  }

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

    data <- data.frame(
      TP = inputs$TP,
      FN = inputs$FN,
      FP = inputs$FP,
      TN = inputs$TN
    )

  data <- correction(data, value = cc.value, type= cc.type)
  data <- logit.data(data)

  return(data)

}


## continuity correction ----
correction <- function(
  data, value = 0.5,
  type = c("single", "all")
){

  type <- match.arg(type)

  if(type == "all"){

    if(any(c(data$TP,data$FN,data$FP,data$TN) == 0)){

      data$TP <- data$TP + value
      data$FN <- data$FN + value
      data$FP <- data$FP + value
      data$TN <- data$TN + value
    }
  }

  if(type == "single"){

    correction = ((((data$TP == 0)|(data$FN == 0))|(data$FP == 0))| (data$TN == 0)) * value

    data$TP <- correction + data$TP
    data$FN <- correction + data$FN
    data$FP <- correction + data$FP
    data$TN <- correction + data$TN

  }

  return(data)

}


## logit-transform data ----
logit.data <- function(data){

  sens <- with(data, TP/(TP+FN))
  spec <- with(data, TN/(TN+FP))

  y1 <- qlogis(sens)
  y2 <- qlogis(spec)

  v1 <- with(data, 1/TP+1/FN)
  v2 <- with(data, 1/TN+1/FP)

  df <- data.frame(
    # sens = sens,
    # spec = spec,
    y1 = y1,
    y2 = y2,
    v1 = v1,
    v2 = v2)

}








