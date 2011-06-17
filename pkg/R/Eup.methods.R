

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
      slope.para <- rbind(inter, slope.para)
    }
   print(t(slope.para))
   cat("\nAdditive Effects Type: ", as.name(x$additive.effects)," \n")
   cat("\nDimension of the Unobserved Factors:", x$used.dim," \n")
   cat("\nNumber of iterations:", x$Nbr.iteration,".")
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