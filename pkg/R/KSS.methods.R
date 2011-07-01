## Methods ========================================================================================

KSS <- function(formula, ...){ UseMethod("KSS")}

print.KSS <- function(x,...){
   cat("Call:\n")
   print(x$call)

   cat("\nCoeff(s) of the Observed Regressor(s) :\n\n")
   slope.para <- x$beta
    if(x$is.intercept){
      Intercept <- matrix(x$mu, 1, 1)
      colnames(Intercept) <- ""
      rownames(Intercept) <- "(Intercept)"
      Para <- rbind(signif(Intercept,digits=3), signif(slope.para,digits=3))
    }
   print(t(Para))
   cat("\nAdditive Effects Type: ", as.name(x$additive.effects)," \n")
   cat("\nDimension of the Unobserved Factors:", x$fAFactMod$used.fdim," \n")
 }

coef.KSS <- function(x,...){
    if(x$is.intercept)
    Intercept        <- x$mu
    else Intercept   <- NULL
    
    Slope.Coef       <- x$beta
    
    if(x$additive.effects== "individual"| x$additive.effects== "twoways")
    Add.Ind.Eff      <- x$tau
    else Add.Ind.Eff <- NULL

    if(x$additive.effects== "time"| x$additive.effects== "twoways")
    Add.Tim.Eff      <- x$beta.0
    else Add.Tim.Eff <- NULL
    
    Factors       <- x$fAFactMod$factors
    
    Loadings      <- x$fAFactMod$loadings
    
    Heterogeneity <- x$fAFactMod$fitted.values
    
    Factor.Dim    <- x$fAFactMod$used.fdim
    
    coef.list     <- list(
        Intercept     = Intercept,
        Slope.Coef    = Slope.Coef,
        Add.Ind.Eff   = Add.Ind.Eff, 
        Add.Tim.Eff   = Add.Tim.Eff, 
        Factors       = Factors, 
        Loadings      = Loadings, 
        Heterogeneity = Heterogeneity,
        Factor.Dim    = Factor.Dim)  
        
    return(coef.list)
}


summary.KSS <- function(x,...){
  ## Residuals & R2: =======================================================================
  Res.outpt <- round((summary(as.vector(x$residuals))), digits=2)[-4]
  names(Res.outpt) <- c("Min", "1Q", "Median", "3Q", "Max")
  yy <- sum(diag(crossprod(x$Orig.Y - mean(x$Orig.Y))))
  ee <- sum(diag(crossprod(x$residuals)))
  R2 <- 1 - ee/yy
  
  ## Add-Effect-Type: ======================================================================
  eff              <- matrix(x$additive.effects)
  colnames(eff)    <- ""
  rownames(eff)    <- ""

  ## Inference for Coefficients: ============================================================
  if(x$is.intercept){
    Intercept.se   <- sqrt(x$Intercept.V)
    beta.se        <- sqrt(diag(x$beta.V))
    Intercept.zval <- coef(x)$Intercept  / Intercept.se
    beta.zval      <- coef(x)$Slope.Coef / beta.se
    TAB  <- cbind("Estimate" = c(signif(as.numeric(coef(x)$Intercept),             digits=3),
                                 signif(as.numeric(coef(x)$Slope.Coef),            digits=3)),
                  "StdErr"   = c(signif(Intercept.se, digits=3),
                                 signif(beta.se,      digits=3)),
                  "z.value"  = c(signif(as.numeric(Intercept.zval),                digits=3),
                                 signif(as.numeric(beta.zval),                     digits=3)),
                  "Pr(>z)"   = c(signif(as.numeric(2*pnorm(-abs(Intercept.zval))), digits=3),
                                 signif(as.numeric(2*pnorm(-abs(beta.zval))),      digits=3))
                  )  
    rownames(TAB)  <- c("(Intercept)", x$names[2:length(x$names)])
  }else{
    beta.se        <- sqrt(diag(x$beta.V))
    beta.zval      <- coef(x)$Slope.Coef / beta.se
    TAB  <- cbind("Estimate" = signif(as.numeric(coef(x)$Slope.Coef),              digits=3),
                  "StdErr"   = signif(beta.se,                                     digits=3),
                  "z.value"  = signif(as.numeric(beta.zval),                       digits=3),
                  "Pr(>z)"   = signif(as.numeric(2*pnorm(-abs(beta.zval))),        digits=3)
                  )  
    rownames(TAB)  <- x$names[2:length(x$names)]
  }
  
  ## Result: ==================================================================================
  result        <- list(Res.outpt    = Res.outpt,
                        coefficients = TAB,
                        R2           = R2,
                        KSS.obj      = x)                
  class(result) <- "summary.KSS"
  return(result)
}

print.summary.KSS <- function(x, ...){
  ## Call
  cat("Call:\n")
  print(x$KSS.obj$call)
  ## Residuals:
  cat("\nResiduals:\n")
  print(x$Res.outpt)
  cat("\n")
  ## Beta-Coeffs
  cat("\n Slope-Coefficients:\n")
  printCoefmat(x$coefficients)
  
  cat("\nAdditive Effects Type: ",              as.name(x$KSS.obj$additive.effects)," \n")
  cat("\nDimension of the Unobserved Factors:", x$KSS.obj$fAFactMod$used.fdim)
  cat("\nOptimized Factor Dimension:         ", x$KSS.obj$fAFactMod$optimal.fdim," \n") 
  cat("\nResidual standard error:",             signif(x$KSS.obj$sig2.hat, digits=3), "on", 
                                                signif(x$KSS.obj$degrees.of.freedom,digits=3), "degrees of freedom \n")
  cat("Multiple R-squared:",                    signif(x$R2,digits=3),"\n")
}

plot.summary.KSS <- function(x,...){
  if(x$KSS.obj$additive.effects=="time"|x$KSS.obj$additive.effects=="twoways"){
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
