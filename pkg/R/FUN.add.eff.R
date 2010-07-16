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
      
    mu        <- ifelse(y.in.list$I, Ym2WO - Xm2WO %*% beta.hat, 0) ## mu:        overall mean effect
    tau       <- (Ym2WI - Ym2WO) - (Xm2WI - Xm2WO) %*% beta.hat     ## tau:       individual effects
    tmp       <- (Ym2WT - Ym2WO) - (Xm2WT - Xm2WO) %*% beta.hat     ## see section 3.1, paper KSS-2009: 
    theta.bar <-  qr.solve(fAFactMod.obj$factors, tmp)              ## theta.bar: scores regarding to TiVC
    beta.0    <-  fAFactMod.obj$factors %*% theta.bar               ## beta.0:    functional time effects

    result           <- matrix(c(mu, tau, beta.0), 1, 3))
    colnames(result) <- c("mu (overall.mean)", "tau (indiv.eff)", "beta.0 (time.eff.function)")
    return(result)
  }
