## rm(list=ls())
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
    est$names         <- names         #Names of: dependent variable and regressors
    est$fAFactMod     <- fAFactMod.obj #Elements: fitted.values,factors,loadings,resid.sd2,given.fdim,optimal.fdim,used.fdim
    class(est)        <- "KSS" 
    ##=============================================================================================
    return(est)
  }



## Methods ========================================================================================

KSS <- function(x, ...){ UseMethod("KSS")}

print.KSS <- function(x,...){
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

coef.KSS <- function(x,...){
  coef.list  <- list(beta=x$beta, beta.0=x$beta.0, mu=x$mu, tau=x$tau)
  names(coef.list$beta) <- x$names[2:length(x$names)]
  return(coef.list)
}

summary.KSS <- function(x,...){
  ## Residuals:
  Res.outpt <- round((summary(as.vector(x$residuals))), digits=2)[-4]
  names(Res.outpt) <- c("Min", "1Q", "Median", "3Q", "Max")
  ## Add-Effect-Type:
  eff              <- matrix(x$effect)
  colnames(eff)    <- ""
  rownames(eff)    <- ""
  ## Coefficients:
  se   <- sqrt(diag(x$beta.V))
  zval <- coef(x)$beta / se
  TAB  <- cbind("Estimate" = round(as.numeric(coef(x)$beta), digits=2),
                "StdErr"   = round(se, digits=2),
                "z.value"  = round(as.numeric(zval), digits=2),
                "Pr(>z)"   = round(as.numeric(2*pnorm(-abs(zval))), digits=2))  
  rownames(TAB) <- x$names[2:length(x$names)]

  ## Number of Dimension and Effects-Type:
  TAB2 <- rbind("Used Factor Dimension: "        = x$fAFactMod$used.fdim,
                "Estimated Factor Dimension: "   = x$fAFactMod$optimal.fdim,
                "Additive Effects Type: "        = as.character(x$eff))
  colnames(TAB2) <- c("")
  
  ## Result:
  result        <- list(call         = x$call,
                        eff          = eff,
                        coefficients = TAB,
                        nDim.effType = TAB2,
                        Res.outpt    = Res.outpt,
                        KSS.obj      = x)                
  class(result) <- "summary.KSS"
  result
}

print.summary.KSS <- function(x, ...){
  ## Call
  cat("Call:\n")
  print(x$call)
  ## Residuals:
  cat("\nResiduals:\n")
  print(x$Res.outpt)
  cat("\n")
  ## Number of Dimension and Effect-Type
  print(x$nDim.effType)
  ## Beta-Coeffs
  cat("\n Beta-Coefficients\n")
  printCoefmat(x$coefficients)
}

plot.summary.KSS <- function(x,...){
  if(x$KSS.obj$eff=="time"|x$KSS.obj$eff=="twoways"){
    par(mfrow=c(1,3))
    plot.ts(x$KSS.obj$beta.0, main="beta.0", ylab="",...)
  }else{par(mfrow=c(1,2))}
  matplot(x$KSS.obj$fAFactMod$factors,
          main=paste("EstimatedFactor-Structure\n(Used Dimension=",x$KSS.obj$fAFactMod$used.fdim,")"),
          xlab="Time",ylab="", type="l",...)
  matplot(x$KSS.obj$fAFactMod$fitted.values,
          main=paste("Fitted individual \nTime-Trends"),
          xlab="Time",ylab="", type="l",...)
  par(mfrow=c(1,1))
}
## ## ================================================================================================
## ## TEST: ==========================================================================================
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/R/OptDim.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/R/pca.fit.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/R/fAFactMod.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/R/FUN.Pformula.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/R/FUN.add.eff.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/R/fpca.fit.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/R/FUN.with.trans.R")
## source("/home/dom/Dokumente/Uni/Promotion/myRoutines/Generate_FS.R")


## ## create data for FPCA
## library(pspline)
## T   = 100
## N   =  50
## dim =   6

## ## FS-Structure
## FS.obj   <- sim.FS(T = T, N = N, dim=dim, Factors= "sin", AR =c(0,0), ar.sd = 0.25, plot.opt = FALSE)
## FS.obs   <- FS.obj[[1]]

## ## Regressor 1
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
## I.scl  <-  matrix(rep(70, N*T),T,N)


## ## Additive-Effects:
##    ## individual-effects

## add.ind      <- sample(c(1:100),N)
## add.ind      <- add.ind-mean(add.ind)
## add.ind      <- matrix(rep(add.ind,each=T),T,N)
##    ## time-effects
## add.tim.fun  <-  FS.obj[[3]] %*% as.matrix(colMeans(FS.obj[[4]]))*c(1e18,1e18,1e18,1e18)
## add.tim.fun  <-  matrix(rep(add.tim.fun,N),T,N); #matplot(add.tim.fun)
## add.tim.fun  <-  add.tim.fun - mean(add.tim.fun[,1])

## ## Panel-Model with Intercept, Global time trend-function, and const individual effects:
## Y            <- I.scl + add.tim.fun + add.ind + 5 * X1 - 5 * X2 + FS.obs; #matplot(Y)

## ## ## Cigarets-Data Set: ##################################################
## ## library(plm);# ?Cigar
## ## data(Cigar)
## ## T        <- (1992-1962)
## ## N        <- 46

## ##   ## Dependent-Var
## ## l.sales  <- log(matrix(Cigar$sales, T,N)) # log.Cigaret-Sales per Capita
## ##   ##
## ## cpi      <- matrix(Cigar$cpi, T,N)        # Consumer Price Index
## ##   ## Independent-Var
## ## l.r.ndi    <- log(matrix(Cigar$ndi,   T,N)/cpi) # log.real.Disposable Income per Capita
## ## l.r.price  <- log(matrix(Cigar$price, T,N)/cpi) # log.real.Price per Pack of Cigarets
## ## l.r.pimin  <- log(matrix(Cigar$pimin, T,N)/cpi) # log.real.Minimum-Price per Pack of Cigarets in Neighbouring States
## ## ##
## ## par(mfrow=c(1,4))
## ## matplot(l.sales,main="l.Sales"); 
## ## matplot(l.r.ndi,main="l.r.ndi");
## ## matplot(l.r.price,main="l.r.Price");
## ## matplot(l.r.pimin,main="l.r.pimin");
## ## par(mfrow=c(1,1))

## ## ## #########################################################################


## ## ## Effecs: "None with Intercept" 
## none.intc.obj      <- KSS(formula=Y       ~ X1          + X2,
##                           effect = "none", dim.crit = "KSS.C1"); #str(none.intc.obj)
## ## Cigs.none.intc.obj <- KSS(formula=l.sales ~ l.r.ndi + l.r.price + l.r.pimin,
## ##                           effect = "none", dim.crit = "KSS.C1", factor.dim=2); #str(Cigs.none.intc.obj)
## summary(none.intc.obj); plot(summary(none.intc.obj))
## ## summary(Cigs.none.intc.obj); # plot(summary(Cigs.none.intc.obj))
## ## check-plot:
## par(mfrow=c(1,3))
## matplot(Y)
## matplot(none.intc.obj$fitted.values)
## matplot(none.intc.obj$residuals)
## par(mfrow=c(1,1))



## ## dieses model ist genau das gleiche wie im KSS-Paper:
## ## Transformation nur mit TimeVaryingConstants "TiVC":
## time.obj      <- KSS(formula=Y       ~-1 + X1      + X2,
##                      effect = "time", dim.crit = "KSS.C1"); #str(time.obj)
## Cigs.time.obj <- KSS(formula=l.sales ~-1 + l.r.ndi + l.r.price + l.r.pimin,
##                      effect = "time", dim.crit = "KSS.C1", factor.dim=5); #str(Cigs.time.obj)

## summary(time.obj); plot(summary(time.obj))
## summary(Cigs.time.obj); plot(summary(Cigs.time.obj))

## ## sollte den Intercept-Parameter ergeben:
## mean(time.obj$beta.0)
## ## check-plot:
## par(mfrow=c(1,3))
## matplot(time.obj$Orig.Y)
## matplot(time.obj$fitted.values)
## matlines(time.obj$beta.0[,1], col="red", lwd=3)
## matplot(time.obj$residuals)
## par(mfrow=c(1,1))

## ## model mit time effekten u individuellen effekten:
## tway.obj      <- KSS(formula=Y       ~-1 + X1 + X2,
##                      effect = "twoways",    dim.crit = "KSS.C1"); # str(tway.obj)
## ## Cigs.tway.obj <- KSS(formula=l.sales ~-1 + l.r.ndi + l.r.price + l.r.pimin,
## ##                      effect = "twoways", dim.crit = "KSS.C1"); # str(Cigs.tway.obj)
## summary(tway.obj)
## ## summary(Cigs.tway.obj); # plot(summary(Cigs.tway.obj))

## round(mean(Cigs.tway.obj$beta.0),digits=2)
## ## check-plot:
## par(mfrow=c(1,2))
## matplot(Y)
## matlines(add.tim.fun, col="red", lwd=3)
## matplot(tway.obj$fitted.values)
## matlines(tway.obj$beta.0[,1], col="red", lwd=3)
## par(mfrow=c(1,1))

## ## model mit time effekten, individuellen effekten und Intercept:
## tway.intcpt.obj <- KSS(formula=Y ~ X1+X2, effect = "twoways",  dim.crit = "KSS.C1"); ## str(tway.intcpt.obj)
## summary(tway.intcpt.obj)
## plot(summary(tway.intcpt.obj))

## tway.intcpt.obj$mu
## mean(tway.intcpt.obj$beta.0)
## ## Intercept:
## I.scl[1,1]
## tway.intcpt.obj$mu
## ## Individual-Effects:
## add.ind[1,]
## mean(tway.intcpt.obj$tau)
## tway.intcpt.obj$tau

## ## check-plot:
## par(mfrow=c(1,4))
## matplot(Y)
## matlines(add.tim.fun, col="red", lwd=3)
## matplot(tway.intcpt.obj$fitted.values)
## matlines(tway.intcpt.obj$beta.0[,1], col="red", lwd=3)
## matplot(tway.intcpt.obj$residuals)
## plot.ts(c(add.tim.fun[,1]-tway.intcpt.obj$beta.0[,1]))
## par(mfrow=c(1,1))




## ## model mit nur mit individuellen effekten:
## indv.obj <- KSS(formula=Y ~ -1+ X1+X2, effect = "individual",  dim.crit = "KSS.C1"); ## str(tway.intcpt.obj)
## Cigs.indv.obj <- KSS(formula=l.sales ~ l.r.ndi + l.r.price + l.r.pimin,
##                      effect = "individual", dim.crit = "KSS.C1", factor.dim=21); ## str(tway.intcpt.obj)
## summary(indv.obj); plot(summary(indv.obj))

## summary(Cigs.indv.obj); plot(summary(Cigs.indv.obj))
## indv.obj$mu
## ## check individual effects:
## indv.obj$tau-70
## add.ind[1,] 
## mean(indv.obj$beta.0)
## ## check-plot:
## par(mfrow=c(1,3))
## matplot(Y)
## matlines(add.tim.fun, col="red", lwd=3)
## matplot(indv.obj$fitted.values)
## ## matlines(indv.obj$beta.0[,1], col="red", lwd=3)
## matplot(indv.obj$residuals)
## par(mfrow=c(1,1))



