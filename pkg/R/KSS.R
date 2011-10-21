## 
KSS.default <- function(formula, consult.dim.crit = FALSE,
                        additive.effects        = c("none", "individual", "time", "twoways"),
                        level = 0.01,
                        factor.dim=NULL,
                        d.max=NULL,
                        sig2.hat = NULL,
                        restrict.mode= c("restrict.factors","restrict.loadings"),...)
  {
    ##===================================================================================
    if(!class(formula)=="formula"){
      stop("\n Argument >>formula<< needs a formula-object like y~x1+... where the elements are matrices.")
    }
    if(!any(additive.effects==c("none", "individual", "time", "twoways"))){
      stop("\n Argument >>effect<< must be one of: \n none, individual, time, twoways")
    }
    if(!is.numeric(level)){
      stop("\n Argument >>alpha<< has to be numeric.")
    }
    ##====================================================================================
    
    ## check "effect" and "dim.crit"
    effect        <- match.arg(additive.effects)
    dim.criterion <- c("PC1", "PC2", "PC3", "IC1", "IC2" , "IC3",
                        "IPC1", "IPC2", "IPC3" , "KSS.C1", "KSS.C2", "ED", "ER", "GR")
    
    ## extract data from formula
    names  <- names(model.frame(formula))
    PF.obj <- FUN.Pformula(formula = formula, effect = effect)
    
    N <- ncol(PF.obj[[1]]$ODM)
    T <- nrow(PF.obj[[1]]$ODM)
    P <- length(PF.obj)-1
    dat.dim 	  <- c(T, N, P)
    is.intercept  <- PF.obj[[1]]$I

    ## *OR*iginal *dat*a
    ORdat    <- sapply(1:(P+1), function(i) PF.obj[[i]]$ODM)
    ## *TR*ansformed *dat*a
    TRdat         <- sapply(1:(P+1), function(i) PF.obj[[i]]$TDM)
    dat.matrix    <- sapply(1:(P+1), function(i) PF.obj[[i]]$TDV) 

    Or.Y     <- ORdat[, 1,       drop = FALSE]					# (TN x 1)
    Or.X     <- ORdat[, 2:(P+1), drop = FALSE]				        # (TN x P)

    TR.Y     <- TRdat[, 1,       drop = FALSE]					# (TN x 1)
    TR.X     <- TRdat[, 2:(P+1), drop = FALSE]				        # (TN x P)
    
    TR.Y.mat <- matrix(TR.Y, T,     N)					        # (T x N)
    TR.X.mat <- matrix(TR.X, T, (N*P))						# (T x NP)

    
    ## smooth.splines with undersmoothing
    ## undersmoothing: 0.8 * GCV-value

    spar.low       <- smooth.Pspline(x = seq.int(1,T), y = TR.Y.mat, spar=0.01, method = 4       )$spar * 0.99
    
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
    ## common-Slope.Coefficients:
    com.slops.0      <- solve(bloc1)%*%bloc2				       # (Px1)

    ## calculate first step residuals and estimate dimension of factor-structure
    Residu.mat    <- matrix((TR.Y - (TR.X %*% com.slops.0)), T, N)

    ## functional pca
    fpca.fit.obj  <- fpca.fit(Residu.mat)

    ## Estimation of Dimension
    Opt.dim.Output <- as.matrix(sapply(dim.criterion, function(dim.criterion){
                                  EstDim(dim.criterion, Obj=Residu.mat, d.max=d.max, factor.dim=factor.dim,
                                         sig2.hat=sig2.hat, level=level)[2]}))
    Opt.dim.Output.Bai <- c(as.numeric(Opt.dim.Output[1:9,1]))
    names(Opt.dim.Output.Bai) <- c("PC1","PC2","PC3","IC1","IC2","IC3","IPC3","IPC2","IPC3")
    Opt.dim.Output.KSS <- c(as.numeric(Opt.dim.Output[10:11,1]))
    names(Opt.dim.Output.KSS) <- c(" KSS.C1","KSS.C1")
    Opt.dim.Output.Onatski <- c(as.numeric(Opt.dim.Output[12,1]))
    names(Opt.dim.Output.Onatski) <- c(" ED")
    Opt.dim.Output.RH <- c(as.numeric(Opt.dim.Output[13:14,1]))
    names(Opt.dim.Output.RH) <- c(" ER","GR")
    if(is.null(factor.dim)&& consult.dim.crit){
      cat("-----------------------------------------------------------\n")
      cat("Results of Dimension-Estimation");cat("\n\n-Bai:\n")
      print(Opt.dim.Output.Bai, quote = FALSE, na.print="");    cat("\n-KSS:\n")
      print(Opt.dim.Output.KSS, quote = FALSE, na.print="");    cat("\n-Onatski:\n")
      print(Opt.dim.Output.Onatski, quote = FALSE, na.print="");cat("\n-RH:\n")
      print(Opt.dim.Output.RH, quote = FALSE, na.print="");     cat("\n")
      cat("-----------------------------------------------------------\n")
      cat("Please, Enter a Dimension:", "\n")
      used.dim <- scan(n=1)
      cat("Used Dimension of unobs. Factor-Structure is:\n", used.dim, "\n")
      cat("-----------------------------------------------------------\n")
    }else {
	if(!is.null(factor.dim)) used.dim <- factor.dim
	else used.dim <- c(as.numeric(Opt.dim.Output[10,1]))
	}
    
    if(used.dim>=1){
      factors       <- fpca.fit.obj$factors[,  1:used.dim, drop= FALSE]
      loadings      <- fpca.fit.obj$loadings[, 1:used.dim, drop= FALSE]
  
      ##==========================================================================
   
      ## re-estimate beta=========================================================
      factor.stract <- tcrossprod(factors, loadings)
      NEW.TR.Y.mat  <- TR.Y.mat - factor.stract
      NEW.TR.Y      <- as.vector(NEW.TR.Y.mat)
      beta          <- qr.solve(TR.X, NEW.TR.Y)
      beta          <- matrix(beta, P,1)
      ##===========================================================================================
    }else{
      factors       <- NULL
      loadings      <- NULL
      beta          <- com.slops.0
      factor.stract <- matrix(0,T,N)
    }
    ## Built up the return object *est*
    
    FUN.add.eff.obj <- FUN.add.eff(PF.obj        = PF.obj,
                                   fpca.fit.obj  = fpca.fit.obj,
                                   beta.hat      = beta)
    
    ## re-estimation of sig2.hat (Paper KSS Section 3.4)==========================================
    YOVc            <- PF.obj[[1]]$TRm$OVc
    XOVc            <- sapply(2:(P+1), function(i)PF.obj[[i]]$TRm$OVc)
    Or.Y_minus_YOVc <- Or.Y - YOVc
    Or.X_minus_XOVc <- Or.X - matrix(rep(XOVc, each=(N*T)), N*T, P)
    
    sig2.hat        <- 1/((N-1)*T) * sum((Or.Y_minus_YOVc - Or.X_minus_XOVc %*% beta - matrix(factor.stract, N*T, 1))^2)

    ## degrees.of.freedom =======================================================================
    degrees.of.freedom <- (T*N - (T+N)*used.dim - P - 
                           is.intercept -
                           N*(effect == "individual"| effect == "twoways") - 
                           T*(effect == "time"      | effect == "twoways"))
    
    ## estimation of Variances  =================================================================
    ## estimation of beta-variance beta.V 
    beta.V <- sig2.hat * solve(bloc1) %*% (t.TR.X.TR.X + t.TR.X.TR.X.smth2 - 2*t.TR.X.TR.X.smth) %*% solve(bloc1)
    ## estimation of Intercept-variance
    
    if(is.intercept){
      Intercept.V <- (sig2.hat + matrix(colMeans(Or.X),1,P) %*% beta.V %*% t(matrix(colMeans(Or.X),1,P)))/(N*T)
    }else{Intercept.V <- NULL}
    
    ## Fitted Values =============================================================================
    fitted.values      <- matrix(c(rep(FUN.add.eff.obj$mu, T*N) + rep(FUN.add.eff.obj$beta.0, N) +
                                   rep(FUN.add.eff.obj$tau, each=T) + Or.X %*% beta)
                                 , T, N) + factor.stract

    ## Return ====================================================================================
    est                    <- vector("list")
    est$dat.matrix         <- dat.matrix
    est$dat.dim            <- dat.dim
    est$slope.para         <- beta
    est$beta.V             <- beta.V
    est$names              <- names         #Names of: dependent variable and regressors
    est$is.intercept       <- is.intercept  #Intercept: TRUE or FALSE
    est$additive.effects   <- effect        #Additive-Effect-Type
    est$Intercept          <- FUN.add.eff.obj$mu     # Intercept
    est$Intercept.V        <- Intercept.V
    est$Add.Ind.Eff        <- FUN.add.eff.obj$tau    # Add. indiv. Effects
    est$Add.Tim.Eff        <- FUN.add.eff.obj$beta.0 # Add. time Effects (beta.0-function)
    est$unob.factors       <- factors
    est$ind.loadings       <- loadings
    est$unob.fact.stru     <- factor.stract
    est$used.dim           <- used.dim
    est$optimal.dim        <- Opt.dim.Output
    est$fitted.values      <- fitted.values
    est$orig.Y             <- matrix(Or.Y, T, N)
    est$residuals          <- est$orig.Y - est$fitted.values
    est$sig2.hat           <- sig2.hat
    est$degrees.of.freedom <- degrees.of.freedom
    est$call               <- match.call()
    class(est)             <- "KSS" 
    ##=============================================================================================
    return(est)
  }
