FUN.add.eff <- function(PF.obj, fAFactMod.obj, beta.hat)
  {
    P         <- length(PF.obj)-1
    y.in.list <- PF.obj[[1]]

# twoways
    
    if(y.in.list$Tr=="twoways"){
      
    #=========================================================================================================
      Ym2WI <- y.in.list$TRm$individual                               # *Y**m*ean *2**W*ay *I*ndividual
      Ym2WT <- y.in.list$TRm$time                                     # *Y**m*ean *2**W*ay *T*ime
      Ym2WO <- y.in.list$OVm                                          # *Y**m*ean *2**W*ay *O*verall           
      Xm2WI <- sapply(2:(P+1), function(i)PF.obj[[i]]$TRm$individual) # *X**m*ean *2**W*ay *I*ndividual
      Xm2WT <- sapply(2:(P+1), function(i)PF.obj[[i]]$TRm$time)       # *X**m*ean *2**W*ay *T*ime
      Xm2WO <- sapply(2:(P+1), function(i)PF.obj[[i]]$OVm)            # *X**m*ean *2**W*ay *O*verall
      
      # Note:
      # "Ym2WO" and "Xm2WO" are garantied to be zero if there is no intercept (see FUN.with.trans()-function)
    #=========================================================================================================
      
      mu        <- ifelse(y.in.list$I, Ym2WO - Xm2WO %*% beta.hat, 0)
      tau       <- (Ym2WI - Ym2WO) - (Xm2WI - Xm2WO) %*% beta.hat
      tmp       <- (Ym2WT - Ym2WO) - (Xm2WT - Xm2WO) %*% beta.hat
      theta.bar <-  qr.solve(fAFactMod.obj$factors, tmp)
      beta.0    <-  fAFactMod.obj$factors %*% theta.bar

      if(y.in.list$I){
        AE        <- structure(list("mu"      = mu,
                                    "beta.0"  = beta.0,
                                    "tau"     = tau),      class = "TwoWwI")
      }else{
        AE        <- structure(list("beta.0"  = beta.0,
                                    "tau"     = tau),      class = "TwoWwoI")
      }
      
# individual or time or none
      
    }else{
      
    #============================================================================================================
      Ym        <- y.in.list$TRm                               # *Y**m*ean (if "time": Tx1, if "individual": Nx1)
      YO        <- y.in.list$OVm                               # *Y**O*verall-mean
      Xm        <- sapply(2:(P+1), function(i)PF.obj[[i]]$TRm) # *X**m*ean (if "time": TxP, if "individual": NxP)
      XO        <- sapply(2:(P+1), function(i)PF.obj[[i]]$OVm) # *X**O*verall-mean
    #============================================================================================================

      mu        <- ifelse(y.in.list$I, YO - XO %*% beta.hat, 0)

      # individual
      if(y.in.list$Tr=="individual"){
        tau <- (Ym - YO) - (Xm - XO) %*% beta.hat
        if(y.in.list$I){
          AE        <- structure(list("mu"  = mu,
                                      "tau" = tau),        class = "IwI")
        }else{
          AE        <- structure(list("tau" = tau),        class = "IwoI")
        } 
      }
      # time
      if(y.in.list$Tr=="time"){
        tmp       <- (Ym - YO) - (Xm - XO) %*% beta.hat
        theta.bar <-  qr.solve(fAFactMod.obj$factors, tmp)
        beta.0    <-  fAFactMod.obj$factors %*% theta.bar
        if(y.in.list$I){
          AE        <- structure(list("mu"     = mu,
                                      "beta.0" = beta.0),  class = "TwI")
        }else{
          AE        <- structure(list("beta.0" = beta.0),  class = "TwoI")
        }
      }
      # none
      if(y.in.list$Tr=="none"){
        if(y.in.list$I){
          AE        <- structure(list("mu"     = mu),      class = "NwI")
        }
      }
    }
    return(AE)
  }
