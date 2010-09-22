#########################################################################################
## FUN.add.eff
## Aim: Calculates Additive effects and Intercept(function)
##
## Takes:
## 1)PF.obj        (return-object of FUN.Pformula(): Original data and Transformed Data)
## 2)fAFactMod.obj (return-object of fAFactMod())
## 3)beta.hat      ((re-)estimated beta coeficients)
##
## Gives a list with:
## mu     (Overallmean)
## tau    (individal effects)
## beta.0 (intercept-function)
##########################################################################################

FUN.add.eff <- function(PF.obj, fAFactMod.obj, beta.hat)
  {
    P         <- length(PF.obj)-1
    y.in.list <- PF.obj[[1]]
      
    ##=========================================================================================================
    YInC  <- y.in.list$TRm$InC                                 ## *Y**In*dividual *C*onstants
    YTiVC <- y.in.list$TRm$TiVC                                ## *Y**Ti*me V*arying *C*onstants
    YOVc  <- y.in.list$TRm$OVc                                 ## *Y**OV*erall *c*onstant          
    XInC  <- sapply(2:(P+1), function(i)PF.obj[[i]]$TRm$InC)   ## *X**In*dividual *C*onstants
    XTiVC <- sapply(2:(P+1), function(i)PF.obj[[i]]$TRm$TiVC)  ## *X**Ti*me V*arying *C*onstants
    XOVc  <- sapply(2:(P+1), function(i)PF.obj[[i]]$TRm$OVc)   ## *X**OV*erall *c*onstant      
    ##=========================================================================================================

    
    mu           <- ifelse(y.in.list$I, YOVc - XOVc %*% beta.hat, 0) ## mu:        overall mean effect
    if(y.in.list$Tr=="individual"|y.in.list$Tr=="twoway"){
      tau        <- c((YInC  - YOVc) - (XInC  - XOVc) %*% beta.hat)  ## tau:       individual effects
    }else{tau    <- 0}    
    if(y.in.list$Tr=="time"|y.in.list$Tr=="twoway"){
      tmp        <- (YTiVC - YOVc) - (XTiVC - XOVc) %*% beta.hat     ## see section 3.1, paper KSS-2009:
      theta.bar  <-  qr.solve(fAFactMod.obj$factors, tmp)            ## theta.bar: scores regarding to TiVC
      beta.0     <-  fAFactMod.obj$factors %*% theta.bar             ## beta.0:    functional time effects
    }else{beta.0 <- 0}
    result    <- list(mu = mu, tau = tau, beta.0 = beta.0)
    return(result)
  }
