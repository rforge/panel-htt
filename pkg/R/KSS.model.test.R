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

KSS.model.test.default <- function(formula,
                                   g.function,
                                   alpha      = 0.05,
                                   ...)# ...  = sign.vec= NULL 
  {
    ## ==================================================================================
    ## effect-only "time"
    effect  <- "time"
    
    ##===================================================================================
    ## check argument formula
    if(!class(formula)=="formula"){
      stop("\n Argument >>formula<< needs a formula-object like y~x1+...")
    }
    if(!is.numeric(alpha)){
      stop("\n Argument >>alpha<< has to be numeric.")
    }
    
    ##====================================================================================
    ## extract data from formula
    names        <- names(model.frame(formula))
    PF.obj       <- FUN.Pformula(formula = formula, effect = effect)
    
    is.intercept <- PF.obj[[1]]$I         # TRUE, if with intercept, FALSE, if not
                                          # (a formula-object with intercept gives: is.intercept==TRUE; see FUN.Pformula)
    N            <- ncol(PF.obj[[1]]$ODM) # Individual-Dimension
    T            <- nrow(PF.obj[[1]]$ODM) # Time-Dimension
    P            <- length(PF.obj)-1      # Number of Regressors (without intercept-term)

    ##====================================================================================
    ## check and evaluate argument >>g.function<<
    if(missing(g.function)){
        ## list-obj with L=1 components
        g.fun   <-  matrix(rep(1, T), T, 1)
        L       <- 1
      }else{
        if(!is.matrix(g.function) & !is.numeric(g.function) & !is.function(g.function) & !is.list(g.function) &
           !is.null(g.function)){
          stop("Wrong argument >>g.function<<. It must be either >>NULL<< (default), a function, a list of functions, a matrix/vector of function-values.")
        }else{
          if(is.matrix(g.function)){# if g.function is a matrix-object
            L     <- dim(g.function)[2]
            g.fun <- g.function           
          }
          if(is.numeric(g.function)&is.null(dim(g.function))){# if g.function is a vector-object
            L <- 1
            g.fun <- g.function                
          }    
          if(is.function(g.function)){# if g.function is a function-object
            g.fun <- try(g.function(c(1:T)), silent=TRUE)
            L     <- 1
            if(class(g.fun)=="try-error"){# check evaluation
              stop("Wrong argument >>g.function<<.\n The given g.function cannot be evaluated at the grid: 1,2,...,T.")
            }
            if(sum(g.fun^2)!=T){# check the size
              stop("Wrong argument >>g.function<<.\n The given g.function has to be of euclidean norm equal to T. Such that: sum(g.function(c(1:T))^2)==T")
            }
            g.fun <- matrix(g.fun, T, L)
          }
          if(is.list(g.function)){
            L     <- length(g.function)
            g.fun <- matrix(NA,    T, L)
            for(l in 1:L){
              if(!is.function(g.function[[l]])){stop("Wrong argument >>g.function<<:\n The function in list-component number ",l," of g.function is not a R function-object.")}
              g.fun[,l] <- tmp <- try(g.function[[l]](c(1:T)), silent=TRUE)
              if(class(tmp)=="try-error"){# check evaluation
                stop(paste("Wrong argument >>g.function<<:\n The function in list-component number ",l," given to the argument g.function cannot be evaluated at the grid: 1,2,...,T."))
              }            
            }
          }
          if(is.null(g.function)){
            warning("No argument >>g.function<< is given. Default factor: g(t)=1 is used.\n This is equivalent to individual (time invariant) effects.")
            ## list-obj with L=1 components
            g.fun   <-  matrix(rep(1, T), T, 1)
            L       <- 1
          }
        }
      }
    ## ==================================================================================
    ## design-matrix and -list (each column of the design.mat (NT x NL) is one list-component)
    design.mat  <- kronecker(diag(N), g.fun)                                    # (NT x NL)
    design.list <- lapply(1:(N*L), function(z, i) matrix(c(z[,i]),T,N), z = design.mat)

    ## ==================================================================================
    ## transform design.list
    g.fun.with.trans.obj  <- sapply(design.list,
                                    function(z) FUN.with.trans(z,
                                                               N            = N,
                                                               T            = T,
                                                               is.intercept = is.intercept,
                                                               effect       = effect), simplify=FALSE)

    
    ## ==================================================================================    
    ## Data Preparations: 
    OR.dat     <- sapply(1:(P+1), function(i) PF.obj[[i]]$ODM)                           # (TN x (P+1))
    ## 
    OR.des     <- sapply(1:(N*L), function(i) g.fun.with.trans.obj[[i]]$ODM)             # (TN x NL)
    
    OR.Y       <- OR.dat[, 1,       drop = FALSE]					        # (TN x 1)
    OR.X       <- OR.dat[, 2:(P+1), drop = FALSE]				                # (TN x P)

    ## ==================================================================================
    ## Orthogonal-Projection Matrix for FWL-Theorem-Application
    M         <- diag(T*N) - (OR.des %*% solve(t(OR.des) %*% OR.des) %*% t(OR.des))         # (TN x TN)

    ## FWL-TRansformations                                                    
    M.OR.Y                 <- M %*% OR.Y                                                # (TN x 1)
    M.OR.X                 <- M %*% OR.X                                                # (TN x P)

    M.OR.Y.mat             <- matrix(M.OR.Y, T, N)                                      # (T x N)
    M.OR.X.mat             <- matrix(M.OR.X, T, N*P)                                    # (T x NP)

    ## ==================================================================================
    ## KSS-Routine

    ## Data
    data.all.mat  <- cbind(M.OR.Y.mat,M.OR.X.mat)

    ## Data-Transformation:
    model.in.list <- lapply(1:(P+1), function(z, i) z[,seq((i-1)*N+1,i*N)], z = data.all.mat)
    data.in.list  <- sapply(model.in.list,
                            function(z) FUN.with.trans(z,
                                                       N            = N,
                                                       T            = T,
                                                       is.intercept = is.intercept,
                                                       effect       = effect),
                            simplify=FALSE)

    M.TR.dat     <- sapply(1:(P+1), function(i) data.in.list[[i]]$TDM)

    M.TR.Y       <- M.TR.dat[,1,         drop = FALSE]
    M.TR.X       <- M.TR.dat[,2:(P+1),   drop = FALSE]
    M.TR.Y.mat   <- matrix(M.TR.Y,T,N)
    M.TR.X.mat   <- matrix(M.TR.X,T,N*P)
    
    spar.low         <- smooth.Pspline(x = seq.int(1,T), y = M.TR.Y.mat, method = 3)$spar * 0.8
    
    M.TR.Y.mat.smth  <- smooth.Pspline(x = seq.int(1,T), y = M.TR.Y.mat,      spar = spar.low)$ysmth       #(T x N)    
    M.TR.X.mat.smth  <- smooth.Pspline(x = seq.int(1,T), y = M.TR.X.mat,      spar = spar.low)$ysmth       #(T x NP)
    M.TR.X.mat.smth2 <- smooth.Pspline(x = seq.int(1,T), y = M.TR.X.mat.smth, spar = spar.low)$ysmth       #(T x NP)

    ## calculate beta coefficents

    M.TR.Y.smth           <- matrix(M.TR.Y.mat.smth,  nrow= (N*T), ncol = 1)	       # (TN x 1)
    M.TR.X.smth           <- matrix(M.TR.X.mat.smth,  nrow= (N*T), ncol = P)	       # (TN x P)
    M.TR.X.smth2          <- matrix(M.TR.X.mat.smth2, nrow= (N*T), ncol = P)	       # (TN x P)

    t.M.TR.X.M.TR.X       <- crossprod(M.TR.X)                                         # (PxP)
    t.M.TR.X.M.TR.X.smth  <- crossprod(M.TR.X, M.TR.X.smth)		               # (PxP)
    t.M.TR.X.M.TR.X.smth2 <- crossprod(M.TR.X, M.TR.X.smth2)		               # (PxP)
    
    t.M.TR.X.M.TR.Y       <- crossprod(M.TR.X, M.TR.Y)     		               # (Px1)
    t.M.TR.X.M.TR.Y.smth  <- crossprod(M.TR.X, M.TR.Y.smth)   		               # (Px1)
    
    M.bloc1               <- t.M.TR.X.M.TR.X - t.M.TR.X.M.TR.X.smth     	       # (PxP)
    M.bloc2               <- t.M.TR.X.M.TR.Y - t.M.TR.X.M.TR.Y.smth                    # (Px1)
    
    M.com.slops.0         <- solve(M.bloc1) %*% M.bloc2				       # (Px1)

    
    ## calculate first step residuals 
    M.Residu.mat    <- matrix((M.TR.Y - M.TR.X %*% M.com.slops.0), T, N)

    ## functional principal component analysis
    M.fpca.fit.obj  <- fpca.fit(dat        = M.Residu.mat,
                             spar.low      = NULL,  # if NULL: 0.8*spar.opt (GCV)
                             given.d       = NULL,  # if NULL: max.rk <- ifelse(neglect.neg.ev,nbr.pos.ev, length(Eval))
                             restrict.mode = "restrict.factors",
                             allow.dual    = TRUE)
    
    ## ===================================================================================
    ## Test Statistik (H_O: dimension of factor structure is == 0)
    
    result  <- KSS.OptDim(Obj      = M.fpca.fit.obj,
                          criteria = c("KSS.C1", "KSS.C2"),
                          sig2.hat = NULL,
                          alpha    = alpha,
                          spar.low = NULL,
                          d.max    = NULL)[[2]]
    result$L       <- L
    if(missing(g.function)){
      result$print <- "TCIV"
    }else{
      if(is.null(g.function)){
        result$print <- "TCIV"
      }else{result$print <- NULL}
    }
    ## ====================================================================================
    class(result)        <- "KSS.model.test" 
    return(result)
  }


## Methods ========================================================================================

KSS.model.test <- function(x,...) UseMethod("KSS.model.test")

print.KSS.model.test <- function(x,...){
  cat("----------------------------------------------\n")
  cat("Test for pre-specified factor model\n")
  cat("----------------------------------------------\n")
  if(is.null(x$print)){
    cat("H0: The given factor model is true.\n\n")
  }else{
    if(x$print=="TCIV"){
      cat("H0: Time-Invariant Indiv.-Effects are TRUE.\n\n")
    }
  }
  outp        <- c(x$Test.Stat, x$p.value, x$crit.value, x$sig.level)
  cat(paste("Number of pre-specified factors: L = ", x$L,"\n\n",sep=""))
  names(outp) <- c("Test-Statistic", "p-value", "crit.-value", "sig.-level")
  print(outp)
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

## #######################################################################################

## ## Model-Test-Checks:

## ## 1) TRUE: Individual Effects:
## Y            <- I.scl + add.ind + 5 * X1 - 5 * X2 + matrix(rnorm(T*N),T,N)
## KSS.model.test(formula=Y ~ X1 + X2, g.function=NULL)
## KSS.model.test(formula=Y ~ X1 + X2)
## ## wrong factros:
## gf <- FS.obj[[3]][,1]
## KSS.model.test(formula=Y ~ X1 + X2, g.function=gf)

## ## 2) TRUE: Individual + Time Effects:
## Y            <- I.scl +add.tim.fun + add.ind + 5 * X1 - 5 * X2 + matrix(rnorm(T*N),T,N)
## KSS.model.test(formula=Y ~ X1 + X2, g.function=NULL)
## KSS.model.test(formula=Y ~ X1 + X2)
## ## wrong factros:
## gf <- FS.obj[[3]][,1]
## KSS.model.test(formula=Y ~ X1 + X2, g.function=gf)




## ## 2) TRUE: 6-dimensional FS:
## Y            <- I.scl + 5 * X1 - 5 * X2 + FS.obs

## g.fun1 <- FS.obj[[3]][,1]
## g.fun2 <- FS.obj[[3]][,2]
## g.fun3 <- FS.obj[[3]][,3]
## g.fun4 <- FS.obj[[3]][,4]
## g.fun5 <- FS.obj[[3]][,5]
## g.fun6 <- FS.obj[[3]][,6]

## gf     <- cbind(g.fun1,g.fun2,g.fun3,g.fun4,g.fun5,g.fun6)
## KSS.model.test(formula=Y ~ X1 + X2, g.function=gf)

## gf     <- cbind(g.fun1,g.fun2,g.fun3,g.fun4,g.fun5)
## KSS.model.test(formula=Y ~ X1 + X2, g.function=gf)

## gf     <- cbind(g.fun1,g.fun2,g.fun3,g.fun4)
## KSS.model.test(formula=Y ~ X1 + X2, g.function=gf)
## ##





