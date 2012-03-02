########################################################################################
## Model-Test (Kneip Sickles Song 2009 Section 3.4)  
##
## Takes:
## formula    = formula-object,
## effect     = c("none", "individual", "time", "twoways")
## g.function = NULL, numeric, matrix, function-obj or a list of function-objects
## alpha      = numeric
##
## Gives:
## Test-Statistic, p.value, crit.value and sig.level
##
## status 12.10.10: not completed
########################################################################################

## rm(list=ls())

ADDvsINT.default <- function(obj,
                             H0.ADD.effects = c("none", "individual", "time", "twoways"),
                             level = 0.05)
  {
    ##===================================================================================
    ## check arguments
    if(!class(obj)=="KSS" & !class(obj)=="Eup"){
      stop("\n Argument >>obj<< needs a either a KSS-object or an Eup-obj.")
    }
    if(!any(H0.ADD.effects==c("none", "individual", "time", "twoways"))){
      stop("\n Argument >>effect<< must be one of: \n none, individual, time, twoways")
    }
    if(!is.numeric(level)){
      stop("\n Argument >>level<< has to be numeric.")
    }
    ## ==================================================================================
    ## effects
    effect       <- match.arg(H0.ADD.effects)
    if(class(obj)=="KSS"){
      ##===================================================================================
      ## calculate residuals
      Residu.mat   <- KSS.default(formula=obj$formula, additive.effects=effect, level=level, factor.dim=0)$residuals
      ## functional principal component analysis
      fpca.fit.obj <- fpca.fit(dat           = Residu.mat,
                               given.d       = NULL,  # if NULL: max.rk <- ifelse(neglect.neg.ev,nbr.pos.ev, length(Eval))
                               restrict.mode = "restrict.factors",
                               allow.dual    = TRUE)   
      ## ===================================================================================
      ## Test Statistik (H_O: dimension of factor structure is == 0)
      
      result  <- KSS.OptDim(Obj      = fpca.fit.obj,
                            criteria = c("KSS.C1"),
                            sig2.hat = NULL,
                            alpha    = level,
                            d.max    = NULL)[[2]]
      result$print <- H0.ADD.effects
    }
    ## ====================================================================================
    class(result)        <- "ADDvsINT" 
    return(result)
  }


## Methods ========================================================================================

ADDvsINT <- function(obj,
                     H0.ADD.effects = c("none","individual", "time", "twoways"),
                     level = 0.05){
  UseMethod("ADDvsINT")
}

print.ADDvsINT <- function(x,...){
  cat("----------------------------------------------\n")
  cat("Testing Additive vs. Interactive Effects\n")
  cat("----------------------------------------------\n")
  cat(paste("H0: There are only additive ",x$print,"-effects\n\n",sep=""))
  outp        <- c(x$Test.Stat, signif(as.numeric(x$p.value),digits=3), x$crit.value, x$sig.level)
  names(outp) <- c("Test-Statistic", "p-value", "crit.-value", "sig.-level")
  print(outp)
}

## ## ================================================================================================
## ## TEST: ==========================================================================================
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/R/OptDim.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/R/pca.fit.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/R/FUN.Pformula.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/R/FUN.add.eff.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/R/fpca.fit.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/R/FUN.with.trans.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/R/KSS.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/R/KSS.methods.R")
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

## ## #######################################################################################

## ## Model-Test-Check:

## ## 1) TRUE: Individual Effects:
## Y            <- I.scl + add.ind + 5 * X1 - 5 * X2 + matrix(rnorm(T*N),T,N)


## ## TEST:
## KSS.obj <- KSS.default(formula=Y ~ X1 + X2, additive.effects="twoways")

## ADDvsINT(obj=KSS.obj, H0.ADD.effects="none")
## ADDvsINT(obj=KSS.obj, H0.ADD.effects="individual")
## ADDvsINT(obj=KSS.obj, H0.ADD.effects="time")
## ADDvsINT(obj=KSS.obj, H0.ADD.effects="twoways")




## load("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel_htt_Arbeitskopie/pkg/data/DGP1.rda")
## N <- 60
## T <- 30

## ## Observed variables:
## Y  <- matrix(DGP1$Y,  nrow=T, ncol=N)
## X1 <- matrix(DGP1$X1, nrow=T, ncol=N)
## X2 <- matrix(DGP1$X2, nrow=T, ncol=N)

## ## Unobserved additive individual-effects:
## add.ind      <- sample(c(1:100),N)
## add.ind      <- add.ind-mean(add.ind)
## add.ind      <- matrix(rep(add.ind,each=T),T,N)

## ## Panel-Data-Generation: 
## Y            <-  add.ind + 5 * X1 - 5 * X2 + matrix(rnorm(T*N),T,N)

## ## KSS-Estimation:
## KSS.obj <- KSS.default(formula=Y ~ X1 + X2)

## #############################################
## ## Test                                    ##
## #############################################

## ## Cannot be rejected
## ADDvsINT(obj=KSS.obj, H0.ADD.effects="individual")
## ## Can be rejected
## ADDvsINT(obj=KSS.obj, H0.ADD.effects="time")

