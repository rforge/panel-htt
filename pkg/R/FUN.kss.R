## rm(list=ls())
FUN.kss.default <- function(formula,
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
                               dim.crit    = dim.crit)
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
    est$sig2.hat <- 1/((N-1)*T) * sum((TR.Y - TR.X %*% beta - matrix(factor.stract, N*T, 1))^2)
    
    ## estimation of beta-variance beta.V ========================================================
    est$beta.V <- est$sig2.hat * solve(bloc1) %*% t.TR.X.TR.X + t.TR.X.TR.X.smth2 - 2*t.TR.X.TR.X.smth %*% solve(bloc1)
    
    ##============================================================================================
    est$fitted.values <- matrix(rep(est$mu, T*N) + rep(est$beta.0, N) + rep(est$tau, each=T) +
                                Or.X %*% beta    + matrix(factor.stract, N*T, 1),
                                T, N)
    est$residuals     <- TR.Y.mat - est$fitted.values
    est$beta          <- beta
    est$call          <- match.call()
    est$effect        <- effect
    est$names         <- names         # names of: dependent variable and regressors
    est$fAFactMod     <- fAFactMod.obj #Elements: fitted.values,factors,loadings,resid.sd2,given.fdim,optimal.fdim,used.fdim
    class(est)        <- "FUN.kss" 
    ##=============================================================================================
    return(est)
  }



## Methods ========================================================================================

FUN.kss <- function(x,...) UseMethod("FUN.kss")

print.FUN.kss <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\nAdditive Effects Type:")
  eff              <- matrix(x$effect)
  colnames(eff)    <- ""
  rownames(eff)    <- ""
  print(eff)
  cat("\nCoeffs of Regressors :\n")
  colnames(x$beta) <- "Beta"
  rownames(x$beta) <- x$names[2:length(x$names)]
  print(x$beta)
}

coef.FUN.kss <- function(x,...){
  coef.list  <- list(beta=x$beta, beta.0=x$beta.0, mu=x$mu, tau=x$tau)
  names(coef.list$beta) <- x$names[2:length(x$names)]
  return(coef.list)
}

summary.FUN.kss <- function(x,...){
  se   <- sqrt(diag(x$beta.V))
  zval <- coef(x)$beta / se
  TAB  <- cbind("Estimate" = as.numeric(coef(x)$beta),
                "StdErr"   = se,
                "z.value"  = as.numeric(zval),
                "Pr(>z)"   = as.numeric(2*pnorm(-abs(zval))))  
  rownames(TAB) <- x$names[2:length(x$names)]
  result        <- list(call=x$call , coefficients=TAB, obj=x)                
  class(result) <- "summary.FUN.kss"
  result
}

print.summary.FUN.kss <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nAdditive Effects Type:")
  eff              <- matrix(x$obj$effect)
  colnames(eff)    <- ""
  rownames(eff)    <- ""
  print(eff)
  cat("Used Dimension:")
  ud <- matrix(x$obj$fAFactMod$used.fdim)
  colnames(ud) <- ""; rownames(ud) <- "";
  print(ud)
  cat("\n")
  printCoefmat(x$coefficients)
}

plot.summary.FUN.kss <- function(x,...){
  if(x$obj$effect=="time"){
    par(mfrow=c(1,3))
    plot.ts(x$obj$beta.0, main="beta.0", ylab="",...)
  }else{par(mfrow=c(1,2))}
  matplot(x$obj$fAFactMod$factors,
          main=paste("EstimatedFactor-Structure\n(Used Dimension=",x$obj$fAFactMod$used.fdim,")"),
          xlab="Time",ylab="", type="l",...)
  matplot(x$obj$fAFactMod$fitted.values,
          main=paste("Fitted individual \nTime-Trends"),
          xlab="Time",ylab="", type="l",...)
  par(mfrow=c(1,1))
}
## ================================================================================================
## TEST: ==========================================================================================
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel-htt/pkg/R/OptDim.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel-htt/pkg/R/pca.fit.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel-htt/pkg/R/fAFactMod.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel-htt/pkg/R/FUN.Pformula.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel-htt/pkg/R/FUN.add.eff.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel-htt/pkg/R/fpca.fit.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel-htt/pkg/R/FUN.with.trans.R")
## source("/home/dom/Dokumente/Uni/Promotion/myRoutines/Generate_FS.R")
## # create data for FPCA
## library(pspline)
## T   = 100
## N   =  50
## dim = 4

## ## FS-Structure
## FS.obj   <- sim.FS(T = T, N = N, dim=dim, Factors= "sin", AR =c(0,0), ar.sd = 0.25)
## FS.obs   <- FS.obj[[1]]

## ## Regressor 2
## X1 <- matrix(NA, T, N)
## for(i in 1:N){
##   X1[,i]       <- seq(1,rnorm(1)*10,length.out=T)+rnorm(T)
## }
## ## Regressor 2
## X2 <- matrix(NA, T, N)
## for(i in 1:N){
##   X2[,i]       <- seq(5,rnorm(1)*1,length.out=T)+rnorm(T,sd=0.75)
## }

## ## Intercept-Scalar
## I.scl  <-  matrix(rep(rnorm(1), N*T),T,N)


## ## Additive-Effects
##  # individual
## add.ind      <- matrix(rep(rnorm(N),each=T),T,N)
##  # time
## add.tim.fun  <-  FS.obj[[3]] %*% as.matrix(colMeans(FS.obj[[4]]))*c(1e20,1e20,1e20,1e20)
## add.tim.fun  <-  matrix(rep(add.tim.fun,N),T,N)

## ## Panel-Model:
## Y        <- add.tim.fun + add.ind + 5 * X1 - 5 * X2 + FS.obs

## ## dieses model ist genau das gleich wie im KSS-Paper:
## ## Transformation nur mit TiVC (i.e.kein overall mean der einzelnen variablen abgezogen):
## time.obj <- FUN.kss(formula=Y~-1+X1+X2, effect = "time",       dim.crit = "KSS.C1"); str(time.obj)

## ## dieses model verzerrt die beta-schätzer stark:
## tway.obj <- FUN.kss(formula=Y~-1+X1+X2, effect = "twoways",    dim.crit = "KSS.C1"); str(tway.obj)

## ## aber man könnte die konstanten individuellen effecte über die scores zurückholen.
## ## gleiches gilt für scalar-intercept.

## none.obj <- FUN.kss(formula=Y~-1+X1+X2, effect = "none",       dim.crit = "KSS.C1"); str(none.obj)
## none.obj <- FUN.kss(formula=Y~   X1+X2, effect = "none",       dim.crit = "KSS.C1"); str(none.obj)

## time.obj <- FUN.kss(formula=Y~-1+X1+X2, effect = "time",       dim.crit = "KSS.C1"); str(time.obj)
## plot(time.obj$beta.0)
## time.obj <- FUN.kss(formula=Y~   X1+X2, effect = "time",       dim.crit = "KSS.C1"); str(time.obj)


## indi.obj <- FUN.kss(formula=Y~   X1+X2, effect = "individual", dim.crit = "KSS.C1"); str(indi.obj)



## none.obj
## indi.obj
## time.obj
## tway.obj

## summary(none.obj)
## summary(indi.obj)
## summary(time.obj)
## summary(tway.obj)

## plot(summary(time.obj))
