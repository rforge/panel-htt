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
    est$fitted.values <- matrix(rep(est$mu, T*N) + rep(est$beta.0, N) + rep(est$tau, each=T) + TR.X %*% beta, T, N)
    est$residuals     <- matrix(TR.Y, T, N) - est$fitted.values
    est$call          <- match.call()
    class(est)        <- "phtt.kss" 
    ##=============================================================================================
    return(est)
  }

