

## Methods ========================================================================================

Eup <- function(formula,
		additive.effects = c("none", "individual",
				     "time", "twoways"),
		dim.criterion	 = c("PC1", "PC2", "PC3", "IC1","IC2" , "IC3", "IPC1", "IPC2","IPC3" , "ED"),
		d.max            = NULL,
		sig2.hat         = NULL,
		factor.dim       = NULL,
		double.iteration = TRUE,
		start.beta       = NULL,
		max.iteration    = 500,
		convergence      = 1e-6,
		restrict.mode    = c("restrict.factors","restrict.loadings"),
                ...){
  UseMethod("Eup")
}

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

summary.Eup <- function(object,...){
  ## Residuals:
  Res.outpt <- round((summary(as.vector(object$residuals))), digits=2)[-4]
  names(Res.outpt) <- c("Min", "1Q", "Median", "3Q", "Max")
  yy <- sum(diag(crossprod(object$orig.Y - mean(object$orig.Y))))
  ee <- sum(diag(crossprod(object$residuals)))
  R2 <- 1 - ee/yy
  #R2a <- 1 - (ee/object$degree.of.freedom)/(yy - 1)
  
  ## Add-Effect-Type:
  eff              <- matrix(object$additive.effects)
  colnames(eff)    <- ""
  rownames(eff)    <- ""
  
  ## Coefficients:
  TAB  <-  Eup.inference(Eup.Obj=object)$inf.result
  
  ## Result:
  result        <- list(Res.outpt    = Res.outpt,
                        coefficients = TAB,
                        R2 = R2,
                        Eup.obj      = object)                
  class(result) <- "summary.Eup"
  result
}


print.summary.Eup <- function(object, ...){
  ## Call
  cat("Call:\n")
  print(object$Eup.obj$call)
  ## Residuals:
  cat("\nResiduals:\n")
  print(object$Res.outpt)
  cat("\n")
  ## Beta-Coeffs
  cat("\n Slope-Coefficients:\n")
  printCoefmat(object$coefficients)
  
  cat("\nAdditive Effects Type: ", as.name(object$Eup.obj$additive.effects)," \n")
  cat("\nDimension of the Unobserved Factors:", object$Eup.obj$used.dim," \n")
  #cat("\nOptimized Factor Dimension:         ", object$Eup.obj$optimal.dim," \n")
  
  cat("\nResidual standard error:", object$Eup.obj$sig2.hat, "on", 
            object$Eup.obj$degree.of.freedom, "degrees of freedom, ", "\nR-squared:", signif(object$R2,digits=3),"\n")
}

plot.summary.Eup <- function(object,...){
  if(is.null(object$Eup.obj$unob.factors) & object$Eup.obj$additive.effects=="none"){
    stop("Neither an estimated factor structure nor additive effects to plot.")
  }
  if(!is.null(object$Eup.obj$unob.factors)){
    if(object$Eup.obj$additive.effects=="none"){
      par(mfrow=c(1,2))
      matplot(object$Eup.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",object$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="l",...)
      matplot(object$Eup.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
    if(object$Eup.obj$additive.effects=="time"){
      par(mfrow=c(1,3))
      plot.ts(object$Eup.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      matplot(object$Eup.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",object$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="l",...)
      matplot(object$Eup.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
    if(object$Eup.obj$additive.effects=="twoways"){
      par(mfrow=c(1,4))
      plot.ts(object$Eup.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      matplot(matrix(rep(object$Eup.obj$Add.Ind.Eff,each=object$Eup.obj$dat.dim[1]),
                     nrow=object$Eup.obj$dat.dim[1],ncol=object$Eup.obj$dat.dim[2]),
                     main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
      matplot(object$Eup.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",object$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="l",...)
      matplot(object$Eup.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
    if(object$Eup.obj$additive.effects=="individual"){
      par(mfrow=c(1,3))
      matplot(matrix(rep(object$Eup.obj$Add.Ind.Eff,each=object$Eup.obj$dat.dim[1]),
                     nrow=object$Eup.obj$dat.dim[1],ncol=object$Eup.obj$dat.dim[2]),
                     main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
      matplot(object$Eup.obj$unob.factors,
            main=paste("Estimated Factors\n(Used Dimension = ",object$Eup.obj$used.dim,")",sep=""),
            xlab="Time",ylab="", type="l",...)
      matplot(object$Eup.obj$unob.fact.stru,
            main=paste("Estimated Factor-Structure"),
            xlab="Time",ylab="", type="l",...)
    par(mfrow=c(1,1))
    }
  }else{
    if(object$Eup.obj$additive.effects=="time"){
      par(mfrow=c(1,1))
      plot.ts(object$Eup.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      par(mfrow=c(1,1))
    }
    if(object$Eup.obj$additive.effects=="twoways"){
      par(mfrow=c(1,2))
      plot.ts(object$Eup.obj$Add.Tim.Eff, main="Additive Time Effect", ylab="",xlab="Time",...)
      matplot(matrix(rep(object$Eup.obj$Add.Ind.Eff,each=object$Eup.obj$dat.dim[1]),
                     nrow=object$Eup.obj$dat.dim[1],ncol=object$Eup.obj$dat.dim[2]),
              main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
      par(mfrow=c(1,1))
    }
    if(object$Eup.obj$additive.effects=="individual"){
      par(mfrow=c(1,1))
      matplot(matrix(rep(object$Eup.obj$Add.Ind.Eff,each=object$Eup.obj$dat.dim[1]),
                     nrow=object$Eup.obj$dat.dim[1],ncol=object$Eup.obj$dat.dim[2]),
                     main="Additive Individual Effects", ylab="",xlab="Time", type="l", ...)
    par(mfrow=c(1,1))
    }
  }
}

