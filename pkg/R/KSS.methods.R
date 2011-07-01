## Methods ========================================================================================

KSS <- function(formula, ...){ UseMethod("KSS")}

print.KSS <- function(x,...){
   cat("Call:\n")
   print(x$call)

   cat("\nCoeff(s) of the Observed Regressor(s) :\n\n")
   slope.para <- x$beta
    if(x$is.intercept){
      inter <- matrix(x$mu, 1, 1)
      colnames(inter) <- ""
      rownames(inter) <- "(Intercept)"
      slope.para <- rbind(signif(inter,digits=3), signif(slope.para,digits=3))
    }
   print(t(slope.para))
   cat("\nAdditive Effects Type: ", as.name(x$effect)," \n")
   cat("\nDimension of the Unobserved Factors:", x$fAFactMod$used.fdim," \n")
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
