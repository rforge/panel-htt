

## Methods ========================================================================================

Eup <- function(formula, ...){ UseMethod("Eup")}

 print.Eup <- function(x,...){
   cat("Call:\n")
   print(x$call)

   cat("\nCoeff(s) of the Observed Regressor(s) :\n\n")
   slope.para <- x$slope.para
    if(x$is.intercept){
      inter <- matrix(x$Intercept, 1, 1)
      colnames(inter) <- ""
      rownames(inter) <- "(Intercept)"
      slope.para <- rbind(signif(inter,digits=3), signif(slope.para,digits=3))
      slope.para <- round(slope.para, 3)
    }
   print(t(slope.para))
   cat("\nAdditive Effects Type: ", as.name(x$additive.effects)," \n")
   cat("\nDimension of the Unobserved Factors:", x$used.dim," \n")
   cat("\nNumber of iterations:", x$Nbr.iteration,".\n")
 }


coef.Eup <- function(x,...){
    if(x$is.intercept)
    Intercept <- x$Intercept
    else Intercept <- NULL
    
    Slope.Coef <- x$slope.para
    
    if(x$additive.effects== "individual"| x$additive.effects== "twoways")
    Add.Ind.Eff <- x$Add.Ind.Eff
    else Add.Ind.Eff <- NULL

    if(x$additive.effects== "time"| x$additive.effects== "twoways")
    Add.Tim.Eff <- x$Add.Tim.Eff
    else Add.Tim.Eff <- NULL
    
    Factors <- x$unob.factors
    
    Loadings <- x$ind.loadings
    
    Heterogeneity <- x$unob.fact.stru
    
    Factor.Dim <- x$used.dim
    
    coef.list <- list(
        Intercept = Intercept,
        Slope.Coef = Slope.Coef,
        Add.Ind.Eff = Add.Ind.Eff, 
        Add.Tim.Eff = Add.Tim.Eff, 
        Factors = Factors, 
        Loadings = Loadings, 
        Heterogeneity = Heterogeneity,
        Factor.Dim = Factor.Dim)  
        
    return(coef.list)
}

residuals.Eup <- function(x, ...){
	x$residuals
	}

summary.Eup <- function(x,...){
  ## Residuals:
  Res.outpt <- round((summary(as.vector(x$residuals))), digits=2)[-4]
  names(Res.outpt) <- c("Min", "1Q", "Median", "3Q", "Max")
  yy <- sum(diag(crossprod(x$orig.Y - mean(x$orig.Y))))
  ee <- sum(diag(crossprod(x$residuals)))
  R2 <- 1 - ee/yy
  #R2a <- 1 - (ee/x$degree.of.freedom)/(yy - 1)
  
  ## Add-Effect-Type:
  eff              <- matrix(x$additive.effects)
  colnames(eff)    <- ""
  rownames(eff)    <- ""
  
  ## Coefficients:
  TAB  <-  Eup.inference(x)$inf.result
  
  ## Result:
  result        <- list(Res.outpt    = Res.outpt,
                        coefficients = TAB,
                        R2 = R2,
                        Eup.obj      = x)                
  class(result) <- "summary.Eup"
  result
}


print.summary.Eup <- function(x, ...){
  ## Call
  cat("Call:\n")
  print(x$Eup.obj$call)
  ## Residuals:
  cat("\nResiduals:\n")
  print(x$Res.outpt)
  cat("\n")
  ## Beta-Coeffs
  cat("\n Slope-Coefficients:\n")
  printCoefmat(x$coefficients)
  
  cat("\nAdditive Effects Type: ", as.name(x$Eup.obj$additive.effects)," \n")
  cat("\nDimension of the Unobserved Factors:", x$Eup.obj$used.dim)
  cat("\nOptimized Factor Dimension:         ", x$Eup.obj$optimal.dim," \n")
  
  cat("\nResidual standard error:", x$Eup.obj$sig2.hat, "on", 
            x$Eup.obj$degree.of.freedom, "degrees of freedom, ", "R-squared:", x$R2,".")
}
plot.summary.Eup <- function(x,...){
  if(is.null(x$Eup.obj$unob.factors) & x$Eup.obj$additive.effects=="none"){
    stop("Neither an estimated factor structure nor additive effects to plot.")
  }
  if(!is.null(x$Eup.obj$unob.factors)){
    if(x$Eup.obj$additive.effects=="none"){
      par(mfrow=c(1,2))
      matplot(x$Eup.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="l",...)
      matplot(x$Eup.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
    if(x$Eup.obj$additive.effects=="time"){
      par(mfrow=c(1,3))
      plot.ts(x$Eup.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      matplot(x$Eup.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="l",...)
      matplot(x$Eup.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
    if(x$Eup.obj$additive.effects=="twoways"){
      par(mfrow=c(1,4))
      plot.ts(x$Eup.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      matplot(matrix(rep(x$Eup.obj$Add.Ind.Eff,each=x$Eup.obj$dat.dim[1]),
                     nrow=x$Eup.obj$dat.dim[1],ncol=x$Eup.obj$dat.dim[2]),
                     main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
      matplot(x$Eup.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="l",...)
      matplot(x$Eup.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
    if(x$Eup.obj$additive.effects=="individual"){
      par(mfrow=c(1,3))
      matplot(matrix(rep(x$Eup.obj$Add.Ind.Eff,each=x$Eup.obj$dat.dim[1]),
                     nrow=x$Eup.obj$dat.dim[1],ncol=x$Eup.obj$dat.dim[2]),
                     main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
      matplot(x$Eup.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",x$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="l",...)
      matplot(x$Eup.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
  }else{
    if(x$Eup.obj$additive.effects=="time"){
      par(mfrow=c(1,1))
      plot.ts(x$Eup.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      par(mfrow=c(1,1))
    }
    if(x$Eup.obj$additive.effects=="twoways"){
      par(mfrow=c(1,2))
      plot.ts(x$Eup.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      matplot(matrix(rep(x$Eup.obj$Add.Ind.Eff,each=x$Eup.obj$dat.dim[1]),
                     nrow=x$Eup.obj$dat.dim[1],ncol=x$Eup.obj$dat.dim[2]),
              main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
      par(mfrow=c(1,1))
    }
    if(x$Eup.obj$additive.effects=="individual"){
      par(mfrow=c(1,1))
      matplot(matrix(rep(x$Eup.obj$Add.Ind.Eff,each=x$Eup.obj$dat.dim[1]),
                     nrow=x$Eup.obj$dat.dim[1],ncol=x$Eup.obj$dat.dim[2]),
                     main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
    par(mfrow=c(1,1))
    }
  }
}

