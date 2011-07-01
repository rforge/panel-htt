## 
KSS.default <- function(formula,
                        effect   = c("none", "individual", "time", "twoways"),
                        dim.crit = c("KSS.C1", "KSS.C2"),
                        alpha    = 0.01,
                        ...)# ...  = sign.vec= NULL 
  {
    ##===================================================================================
    if(!class(formula)=="formula"){
      stop("\n Argument >>formula<< needs a formula-object like y~x1+...")
    }
    if(!any(effect==c("none", "individual", "time", "twoways"))){
      stop("\n Argument >>effect<< must be one of: \n none, individual, time, twoways")
    }
    if(!is.numeric(alpha)){
      stop("\n Argument >>alpha<< has to be numeric.")
    }
    ##====================================================================================
    
    ## check "effect" and "dim.crit"
    effect   <- match.arg(effect)
    dim.crit <- match.arg(dim.crit)
    
    ## extract data from formula
    names  <- names(model.frame(formula))
    PF.obj <- FUN.Pformula(formula = formula, effect = effect)
    
    N <- ncol(PF.obj[[1]]$ODM)
    T <- nrow(PF.obj[[1]]$ODM)
    P <- length(PF.obj)-1
    is.intercept <- PF.obj[[1]]$I

    ## *OR*iginal *dat*a
    ORdat    <- sapply(1:(P+1), function(i) PF.obj[[i]]$ODM)
    ## *TR*ansformed *dat*a
    TRdat    <- sapply(1:(P+1), function(i) PF.obj[[i]]$TDM) 

    Or.Y     <- ORdat[, 1,       drop = FALSE]					# (TN x 1)
    Or.X     <- ORdat[, 2:(P+1), drop = FALSE]				        # (TN x P)

    TR.Y     <- TRdat[, 1,       drop = FALSE]					# (TN x 1)
    TR.X     <- TRdat[, 2:(P+1), drop = FALSE]				        # (TN x P)
    
    TR.Y.mat <- matrix(TR.Y, T,     N)					        # (T x N)
    TR.X.mat <- matrix(TR.X, T, (N*P))						# (T x NP)

    
    ## smooth.splines with undersmoothing
    ## undersmoothing: 0.8 * GCV-value

    spar.low       <- smooth.Pspline(x = seq.int(1,T), y = TR.Y.mat, method = 3       )$spar * 0.8
    
    TR.Y.mat.smth  <- smooth.Pspline(x = seq.int(1,T), y = TR.Y.mat,      spar   = spar.low)$ysmth       #(T x N)    
    TR.X.mat.smth  <- smooth.Pspline(x = seq.int(1,T), y = TR.X.mat,      spar   = spar.low)$ysmth       #(T x NP)
    TR.X.mat.smth2 <- smooth.Pspline(x = seq.int(1,T), y = TR.X.mat.smth, spar   = spar.low)$ysmth       #(T x NP)
    
    ## calculate beta coefficents

    TR.Y.smth        <- matrix(TR.Y.mat.smth,  nrow= (N*T), ncol = 1)	       # (TN x 1)
    TR.X.smth        <- matrix(TR.X.mat.smth,  nrow= (N*T), ncol = P)	       # (TN x P)
    TR.X.smth2       <- matrix(TR.X.mat.smth2, nrow= (N*T), ncol = P)	       # (TN x P)

    t.TR.X.TR.X      <- crossprod(TR.X)       			               # (PxP)
    t.TR.X.TR.X.smth <- crossprod(TR.X, TR.X.smth)		               # (PxP)
    t.TR.X.TR.X.smth2<- crossprod(TR.X, TR.X.smth2)		               # (PxP)
    
    t.TR.X.TR.Y      <- crossprod(TR.X, TR.Y)     		               # (Px1)
    t.TR.X.TR.Y.smth <- crossprod(TR.X, TR.Y.smth)   		               # (Px1)
    
    bloc1            <- t.TR.X.TR.X - t.TR.X.TR.X.smth     		       # (PxP)
    bloc2            <- t.TR.X.TR.Y - t.TR.X.TR.Y.smth     	               # (Px1)
    
    com.slops.0 <- solve(bloc1)%*%bloc2					       # (Px1)


    ## calculate first step residuals and estimate dimension of factor-structure
    Residu.mat    <- matrix((TR.Y - TR.X %*% com.slops.0), T, N)

    fAFactMod.obj <- fAFactMod(dat         = Residu.mat,
                               ## the following args go to FUN.with.trans() and leave Residu.mat unchanged
                               ## add.effects = "none" & demean = FALSE: Because, with.transformations is already done above.
                               demean      = FALSE,
                               add.effects = "none",
                               ## the following args go to KSS.dim.opt() via EstDim()
                               alpha       = alpha,
                               dim.crit    = dim.crit, ...)
    ## *fAFactMod.obj* is a list with: fitted.values,factors,loadings,resid.sd2,given.fdim,optimal.fdim,used.fdim    
    ##==========================================================================


   
    ## re-estimate beta=========================================================
    factor.stract <- tcrossprod(fAFactMod.obj$factors, fAFactMod.obj$loadings)
    NEW.TR.Y.mat  <- TR.Y.mat - factor.stract
    NEW.TR.Y      <- as.vector(NEW.TR.Y.mat)
    beta          <- qr.solve(TR.X, NEW.TR.Y)
    beta          <- matrix(beta, P,1)
    ##===========================================================================================

    ## Built up the return object *est*
    
    est <- FUN.add.eff(PF.obj        = PF.obj,
                       fAFactMod.obj = fAFactMod.obj,
                       beta.hat      = beta)
    
    ## re-estimation of sig2.hat (Paper KSS Section 3.4)==========================================
    YOVc            <- PF.obj[[1]]$TRm$OVc
    XOVc            <- sapply(2:(P+1), function(i)PF.obj[[i]]$TRm$OVc)
    Or.Y_minus_YOVc <- Or.Y - YOVc
    Or.X_minus_XOVc <- Or.X - matrix(rep(XOVc, each=(N*T)), N*T, P)
    
    est$sig2.hat <- 1/((N-1)*T) * sum((Or.Y_minus_YOVc - Or.X_minus_XOVc %*% beta - matrix(factor.stract, N*T, 1))^2)
    
    ## estimation of beta-variance beta.V ========================================================
    est$beta.V <- est$sig2.hat * solve(bloc1) %*% (t.TR.X.TR.X + t.TR.X.TR.X.smth2 - 2*t.TR.X.TR.X.smth) %*% solve(bloc1)
    
    ##============================================================================================
    est$fitted.values <- matrix(rep(est$mu, T*N) + rep(est$beta.0, N) + rep(est$tau, each=T) +
                                Or.X %*% beta    + matrix(fAFactMod.obj$fitted.values,N*T,1),#matrix(factor.stract, N*T, 1),
                                T, N)
    est$Orig.Y        <- matrix(Or.Y, T, N)
    est$residuals     <- est$Orig.Y - est$fitted.values
    est$beta          <- beta
    est$call          <- match.call()
    est$effect        <- effect        #Additive-Effect-Type
    est$is.intercept  <- is.intercept  #Intercept: TRUE or FALSE
    est$names         <- names         #Names of: dependent variable and regressors
    est$fAFactMod     <- fAFactMod.obj #Elements: fitted.values,factors,loadings,resid.sd2,given.fdim,optimal.fdim,used.fdim
    class(est)        <- "KSS" 
    ##=============================================================================================
    return(est)
  }







