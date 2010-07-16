FUN.add.eff <- function(PF.obj, fAFactMod.obj, beta.hat)
  {
    P         <- length(PF.obj)-1
    y.in.list <- PF.obj[[1]]
      
    ##=========================================================================================================
    YInC  <- y.in.list$TRm$InC                                 ## *Y**In*dividual *C*onstants
    YTiVC <- y.in.list$TRm$TiVC                                ## *Y*
    YOVc  <- y.in.list$TRm$OVc                                 ## *Y*          
    XInC  <- sapply(2:(P+1), function(i)PF.obj[[i]]$TRm$InC)   ## *X*
    XTiVC <- sapply(2:(P+1), function(i)PF.obj[[i]]$TRm$TiVC)  ## *X*
    XOVc  <- sapply(2:(P+1), function(i)PF.obj[[i]]$TRm$OVc)   ## *X*      
    ##=========================================================================================================
      
    mu        <- ifelse(y.in.list$I, YOVc - XOVc %*% beta.hat, 0) ## mu:        overall mean effect
    tau       <- (YInC  - YOVc) - (XInC  - XOVc) %*% beta.hat     ## tau:       individual effects
    tmp       <- (YTiVC - YOVc) - (XTiVC - XOVc) %*% beta.hat     ## see section 3.1, paper KSS-2009:
    print("hier")
    theta.bar <-  qr.solve(fAFactMod.obj$factors, tmp)            ## theta.bar: scores regarding to TiVC
    beta.0    <-  fAFactMod.obj$factors %*% theta.bar             ## beta.0:    functional time effects

    result    <- list(mu = mu, tau = tau, beta.0 = beta.0)
    return(result)
  }
