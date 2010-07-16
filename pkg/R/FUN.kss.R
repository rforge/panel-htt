FUN.kss <- function(formula,
                            effect = c("time", "individual", "twoway", "none"),
                            alpha  = 0.05,
                            ...)# ...  = sign.vec= NULL 
  {
    ##===================================================================================
    if(!class(formula)=="formula"){
      stop("\n Argument >>formula<< needs a formula-object like y~x1+...")
    }
    if(!any(effect==c("time", "individual", "twoway", "none"))){
      stop("\n Argument >>effect<< must be one of: \n time, individual, twoway, none")
    }
    if(!is.numeric(alpha)){
      stop("\n Argument >>alpha<< has to be numeric.")
    }
    ##====================================================================================
    
    ## check effect.mode
    effect <- match.arg(effect)
    
    ## extract data from formula
    
    PF.obj <- FUN.Pformula(formula = formula, effect = effect)
    
    N <- ncol(PF.obj[[1]]$ODM)
    T <- nrow(PF.obj[[1]]$ODM)
    P <- length(PF.obj)-1
    
    TRdat    <- sapply(1:(P+1), function(i) PF.obj[[i]]$TDM) 
    
    TR.Y     <- TRdat[, 1,       drop = FALSE]					# (TN x 1)
    TR.X     <- TRdat[, 2:(P+1), drop = FALSE]				        # (TN x P)
    
    TR.Y.mat <- matrix(TR.Y, T,     N)					        # (T x N)
    TR.X.mat <- matrix(TR.X, T, (N*P))						# (T x NP)

    

    
    ## smooth.splines with undersmoothing
    ## undersmoothing: 0.8 * GCV-value

    spar.low       <- smooth.Pspline(x = seq.int(1,T), y = TR.Y.mat, method = 3       )$spar * 0.8
    
    TR.Y.mat.smth  <- smooth.Pspline(x = seq.int(1,T), y = TR.Y.mat, spar   = spar.low)$ysmth       #(T x N)    
    TR.X.mat.smth  <- smooth.Pspline(x = seq.int(1,T), y = TR.X.mat, spar   = spar.low)$ysmth       #(T x NP)
    
    ## calculate beta coefficents

    TR.Y.smth        <- matrix(TR.Y.mat.smth, nrow= (N*T), ncol = 1)	       # (TN x 1)
    TR.X.smth        <- matrix(TR.X.mat.smth, nrow= (N*T), ncol = P)	       # (TN x P)

    t.TR.X.TR.X      <- crossprod(TR.X)       			               # (PxP)
    t.TR.X.TR.X.smth <- crossprod(TR.X, TR.X.smth)		               # (PxP)
    
    t.TR.X.TR.Y      <- crossprod(TR.X, TR.Y)     		               # (Px1)
    t.TR.X.TR.Y.smth <- crossprod(TR.X, TR.Y.smth)   		               # (Px1)
    
    bloc1            <- t.TR.X.TR.X - t.TR.X.TR.X.smth     		       # (PxP)
    bloc2            <- t.TR.X.TR.Y - t.TR.X.TR.Y.smth     	               # (Px1)
    
    com.slops.0 <- solve(bloc1)%*%bloc2					       # (Px1)


    ## calculate first step residuals and estimate dimension of factor-structure
    Residu.mat    <- matrix((TR.Y - TR.X %*% com.slops.0), T, N)

    fAFactMod.obj <- fAFactMod(dat      = Residu.mat,
                               alpha    = alpha,
                               dim.crit = "KSS")
    ##==========================================================================

    ## re-estimate beta=========================================================
    factor.stract <- tcrossprod(fAFactMod.obj$factors, fAFactMod.obj$loadings)  
    NEW.TR.Y.mat  <- TR.Y.mat - factor.stract
    NEW.TR.Y      <- as.vector(NEW.TR.Y.mat)
    beta          <- qr.solve(TR.X, NEW.TR.Y)
    ##==========================================================================

    est <- FUN.add.eff(PF.obj        = PF.obj,
                       fAFactMod.obj = fAFactMod.obj,
                       beta.hat      = beta)
    

    ##=============================================================================================
    est$fitted.values <- rep(est$mu, T*N) + rep(est$beta.0, N) + rep(est$tau, each=T) + TR.X %*% est$beta
    est$residuals     <- TR.Y - fitted.values
    est$call          <- match.call()
    class(est)        <- "phtt.kss" 
    ##=============================================================================================
    return(est)
  }








##  # functional intercept and additive effects      
##     if(effect=="time"){
##     if(PF.obj[[1]]$I){# intercept==T
      
##       #================================================================      
##       x.all.OVm                   <- NULL
##       for(p in 2:(P+1)){x.all.OVm <- c(    x.all.OVm, PF.obj[[p]]$OVm)}
##       x.all.TRm                   <- NULL
##       for(p in 2:(P+1)){x.all.TRm <- cbind(x.all.TRm, PF.obj[[p]]$TRm)}
##       #================================================================
##       mu                          <- PF.obj[[1]]$OVm - x.all.OVm %*% beta
##       if(effect=="time"|effect=="twoway"){
##         tmp                         <- matrix((PF.obj[[1]]$TRm - PF.obj[[1]]$OVm),ncol=1) -
##                                       (x.all.TRm - matrix(rep(x.all.OVm, each=T), nrow=T, ncol=P)) %*% beta
##       }
##       if(effect=="individual"|effect=="twoway"){
##         tau                         <- matrix((PF.obj[[1]]$TRm - PF.obj[[1]]$OVm),ncol=1) -
##                                       (x.all.TRm - matrix(rep(x.all.OVm, each=T), nrow=T, ncol=P)) %*% beta
##       }      
##       theta.bar                   <- qr.solve(fAFactMod.obj$factors, tmp)
##       beta.0                      <- fAFactMod.obj$factors %*% theta.bar
##       # recovering the beta_0-function as in Kneip et al.
##       beta.0                      <- beta.0 + as.numeric(mu)
      
##     }else{# intercept==F
##       #================================================================      
##       x.all.TRm                   <- NULL
##       for(p in 2:(P+1)){x.all.TRm <- cbind(x.all.TRm, PF.obj[[p]]$TRm)}
##       #================================================================      
      
##       tmp     <- matrix(PF.obj[[1]]$TRm, ncol=1) - x.all.TRm %*% beta
      
##       theta.bar <- qr.solve(  fAFactMod.obj$factors, tmp)
##       beta.0    <- fAFactMod.obj$factors %*% theta.bar
      
##     }
##   }








##             raw.f.interc.indv         <- PF.obj[[1]]$OVm + PF.obj[[1]]$TRm$individual - x.all.means %*% beta
##             qr.solve(fAFactMod.obj$factors, 
            
##             f.intercept.indv.smth <- smooth.Pspline(x = basis.t , y = f.intercept.indv, df=df.2, method=2)$ysmth
##             f.intercept.time      <- y.in.list$OVm + y.mean.twoway.time.vec - x.all.mean.twoway.time %*% beta
##             f.intercept.time.smth <- smooth.Pspline(x = basis.t , y = f.intercept.time, df=df.2, method=2)$ysmth
##           }else{
##             f.intercept.indv      <- y.mean.twoway.indv.vec - x.all.mean.twoway.indv %*% beta
##             f.intercept.indv.smth <- smooth.Pspline(x = basis.t , y = f.intercept.indv, df=df.2, method=2)$ysmth
##             f.intercept.time      <- y.mean.twoway.time.vec - x.all.mean.twoway.time %*% beta
##             f.intercept.time.smth <- smooth.Pspline(x = basis.t , y = f.intercept.time, df=df.2, method=2)$ysmth
##           }          
##         }else{
##           if(is.intercept){
##             f.intercept         <- y.in.list$OVm + y.mean.vec - x.all.mean %*% beta
##             f.intercept.smth    <- smooth.Pspline(x = basis.t , y = f.intercept, df=df.2, method=2)$ysmth
##           }else{
##             f.intercept         <- y.mean.vec - x.all.mean %*% beta
##             f.intercept.smth    <- smooth.Pspline(x = basis.t , y = f.intercept, df=df.2, method=2)$ysmth
##           }
##         }
        
## # fitted values

## 	x.beta   <- matrix(x.all.matrix, N*T, P)%*% beta

##         if(y.in.list$T=="twoway"){
##           if(is.intercept){
##             y.fitted <- matrix(y.in.list$OVm, T, N) + matrix(f.intercept.time.smth, T, N) +
##               matrix(f.intercept.indiv.smth, T, N) +  matrix(x.beta, T, N) + factor.stract
##           }else{
##             y.fitted <- matrix(f.intercept.time.smth, T, N) + matrix(f.intercept.indiv.smth, T, N) +
##               matrix(x.beta, T, N) + factor.stract
##           }
##         }else{
##           if(is.intercept){
##             y.fitted <- matrix(y.in.list$OVm, T, N) + matrix(f.intercept.smth, T, N) +
##               matrix(x.beta, T, N) + factor.stract
##           }else{
##             y.fitted <- matrix(f.intercept.smth, T, N) + matrix(x.beta, T, N) + factor.stract
##           }
##         }
                  
##         obj <- list(basis.t, y.matrix, x.all.matrix, N,T,P, df.1, df.2, dim.bf, lambda, eig.func, Residu.mat, Resi.smth.mat)
        
        
## par(mfrow= c(2,3))
## matplot(scores, typ = "l")
## matplot(eig.func, typ = "l")
## matplot(factor.stract, typ= "l")
## plot(f.intercept.smth, typ = "l")
## matplot(y.fitted, typ = "l")
## matplot(y.matrix, typ= "l")
## print(beta)# from simulation: beta = c(-1.5, 5, -1.8)
##         obj


        
