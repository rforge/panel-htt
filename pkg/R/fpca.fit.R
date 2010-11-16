########### spectral variance decomposition and pca fitting if asked ###########################

fsvd.pca.ramsay <- function(Q,
                            allow.dual      = TRUE,
                            given.d         = NULL,
                            calcul.loadings = TRUE,
                            neglect.neg.ev  = FALSE,
                            spar.low        = NULL   # No smoothing is done, if: spar.low = 0
                     ){

  ## extract data information
  nr      <- nrow(Q)
  nc      <- ncol(Q)

  ## save original Q-values 
  Q.non.smth <- Q 

  ## smoothing Q (small degree of undersmoothing)==============================================#
  if(is.null(spar.low)){                                                                       #
    spar.low <- smooth.Pspline(x=seq(0, 1, length.out=nr), y=Q, method = 3       )$spar  * 0.8 #
  }                                                                                            #
  Q          <- smooth.Pspline(x=seq(0, 1, length.out=nr), y=Q, spar   = spar.low)$ysmth       #
  ##===========================================================================================#

  ## For method=Ramsay: dual-matrix or not?
  dual    <- (nr>nc & allow.dual | calcul.loadings )
 
  if(dual){
    Q <- t(Q)
  }

  ## classical approach of Ramsay et al.
  ## trapezoidal rule (for integral-approximation)=============================#
  ## Hier nicht unbedingt notwendig, da                                        #
  ## angenommen wird, dass:                                                    #
  ## -beide indizes i bzw. t bei 1 beginnen!                                   #
  ## -len.Interval==n.discr-1 s.t: (len.Interval)/(n.discr-1) = 1              #
  len.Interval    <- ifelse(dual, nc-1, nr-1)       # dual: N-1, non-dual: T-1 #
  n.discr         <- ifelse(dual, nc,   nr  )       # dual: N  , non-dual: T   #  
  h               <- (len.Interval)/(n.discr-1)                                #
  w               <- c(h/2, rep(h, n.discr-2), h/2)                            #
  cov.mat         <- diag(sqrt(w)) %*% tcrossprod(Q) %*% diag(sqrt(w))         #
  ##===========================================================================#
  
  ## Compute spectral decompostion 
  Spdec           <- eigen(cov.mat, symmetric= TRUE)
  Eval	          <- Spdec[[1]]
  Evec    	  <- Spdec[[2]]

  ## compare rank and given.d  
  nbr.pos.ev      <- length(Eval[Eval > 0])
  max.rk          <- ifelse(neglect.neg.ev, nbr.pos.ev, length(Eval))
  if(is.null(given.d)){
    given.d <- max.rk
  }else{
    if(given.d > max.rk){
      warning(c("The given dimension 'given.d' is larger than the number of positve eigen values."))
    }
    given.d <- min(given.d, max.rk)
  }

  ## compute spectral variance decomposition
  ######################################################
  ## Approximation of the L2-normalization-constraint: #
  ##    ||e.fun_i||_L2 = 1 approximated by:            #
  ##  w*||e.fun_i||_E  = 1                             #
  ######################################################

  ## Left side decomposion
  L                         <- Evec[,0:max.rk , drop= FALSE]   # dual==FALSE: (T x max.rk), dual==TRUE: (N x max.rk)
  ## *fun*ctional approximation
  L.fun                     <- diag(1/(sqrt(w))) %*% L

  ## sqrt of eigenvalues      
  if(!neglect.neg.ev){
    sqr.E                   <- c(sqrt(Eval[Eval > 0]), rep(0, (max.rk - nbr.pos.ev)))
  }else{
    sqr.E	            <- sqrt(Eval[1:max.rk ])
  }

  ## computation of the loadings-parameter (scores)
  if(calcul.loadings){
    S                       <- crossprod(Q, L)[, 0:max.rk , drop= FALSE] # crossprod: t(Q) %*% L
                                                                         # if dual: dim(S)= TxN, if non-dual: dim(S)=NxT
    ## scaling such that: ||R.fun||_E == I
    R                       <- S %*% diag(1/(sqrt(diag(crossprod(S)))))  # if dual: dim(R)= TxN, if non-dual: dim(R)=NxT
    ## approximation to L2-norm ===========================================================================#
    len.Interval.ast    <- ifelse(!dual, nc-1, nr-1)                     # if dual: N-1, if non-dual: T-1  #
    n.discr.ast         <- ifelse(!dual, nc,   nr  )                     # if dual: N  , if non-dual: T    #
    h.ast               <- (len.Interval.ast)/(n.discr.ast-1)                                              #
    w.ast               <- c(h.ast/2, rep(h.ast, n.discr.ast-2), h.ast/2)                                  #
    R.fun               <- diag(1/(sqrt(w.ast))) %*% R                                                     #
    ##=====================================================================================================#

    ## no pca-fitting
    if(((given.d==max.rk)&&!neglect.neg.ev)){
      Q.fit                 <- Q
    ## pca-fitting  
    }else{
      Q.fit                 <- tcrossprod(L.fun[, 0:given.d , drop= FALSE],
                                          L.fun[, 0:given.d , drop= FALSE]) %*% Q
    }
  ## no computation of the loadings-parameter (scores)  
  }else{
    R.fun                   <- NULL
    if(((given.d==max.rk)&&!neglect.neg.ev)){
      Q.fit <- Q
    }else{
      Q.fit                 <- tcrossprod(L.fun[, 0:given.d , drop= FALSE],
                                          L.fun[, 0:given.d , drop= FALSE]) %*% Q
    }
  }

  ## re-convert dimension if dual covariance matrix was used
  if(dual){
    u          <- L.fun
    L.fun      <- R.fun
    R.fun      <- u
    Q.fit      <- t(Q.fit)
    Q          <- t(Q)          # (under)smoothed Q
    Q.non.smth <- t(Q.non.smth) # non-smoothed Q    
  }
    
  ## prepare return-values       
  d.seq <- seq.int(0, (max.rk-1)) # dimension-sequence
  E     <- Eval[1:max.rk]
  sum.e <- sum(E)
  cum.e <- cumsum(E)
  V.d   <- c(sum.e, sum.e-cum.e[-length(cum.e)])

  ## return
  structure(list(L              = L.fun,
                 R              = R.fun,
                 Q.orig         = Q.non.smth,
                 Q.orig.smth    = Q,
                 Q.fit          = Q.fit,
                 spar.low       = spar.low,
                 E              = E,
                 sqr.E          = sqr.E,
                 given.d        = given.d,
                 d.seq          = d.seq,
                 V.d            = V.d,
                 nr             = nr,
                 nc             = nc ,
                 cov.mat        = cov.mat,
                 dual           = dual,
                 fpca.method    = "Ramsay"),
            class   = "fsvd.pca")
}

#######################################################################################
## DUAL-approach (Benko, Härdle, Kneip 2008: "Common functional principal components")#
#######################################################################################
fsvd.pca.kneip <- function(Q,
                           given.d         = NULL,
                           calcul.loadings = TRUE,
                           neglect.neg.ev  = FALSE,
                           spar.low        = NULL   # No smoothing is done, if: spar.low = 0
                           ){
  
  ## extract data information
  nr      <- nrow(Q)
  nc      <- ncol(Q)
  
  M       <- crossprod(Q)                                                             # t(TxN)%*%(TxN)=(NxN)
  ## Weights for non param variance estimator
  ## (Rice 1984 bzw. Hall, Key & Titterington 1990)
  d.0         <-  1/sqrt(1)
  d.1         <- -1/sqrt(1)
  sig2.hat    <- (colSums((d.0*Q[-nr,] + d.1*Q[-1,])^2)) / (nr - 1)                   # (1xN)
  ## Subtract the Variance from the observation-errors
  M           <- M-diag(sig2.hat)                                                     # (NxN)=(NxN)-(NxN)
  
  ## Compute spectral decompostion of the M matrix
  M.Spdec     <- eigen(M, symmetric= TRUE)
  M.Eval      <- M.Spdec[[1]]
  M.Evec      <- M.Spdec[[2]]

  ## Compare rank and given.d  
  nbr.pos.ev      <- length(M.Eval[M.Eval > 0])
  max.rk          <- ifelse(neglect.neg.ev, nbr.pos.ev, length(M.Eval))
  if(is.null(given.d)){
    given.d <- max.rk
  }else{
    if(given.d > max.rk){
      warning(c("The given dimension 'given.d' is larger than the number of positve eigen values."))
    }
    given.d <- min(given.d, max.rk)
  }
  

  R                 <- M.Evec                                                         # (N x max.rk)
 
  ## save original Q-values 
  Q.non.smth <- Q 

  ## smoothing Q (small degree of undersmoothing)================================================#
  if(is.null(spar.low)){                                                                         #
    spar.low <- smooth.Pspline(x=seq(0, 1, length.out=nr), y=Q, method = 3       )$spar  * 0.8   #
  }                                                                                              #
  Q          <- smooth.Pspline(x=seq(0, 1, length.out=nr), y=Q, spar   = spar.low)$ysmth # (TxN) #
  ##=============================================================================================#

  sqrt.inv.eval.mat  <- matrix(rep(1/(sqrt(M.Eval[1:max.rk ])), each=nr), nr, max.rk)
  L                  <- sqrt.inv.eval.mat * Q %*% R                                   # (TxN)%*%(N x max.rk)=(T x max.rk)

  ## no fpca-fitting
    if(((given.d==max.rk)&&!neglect.neg.ev)){
      Q.fit                 <- Q
      ## fpca-fitting  
    }else{
      Q.fit                 <- tcrossprod(L[, 0:given.d , drop= FALSE],
                                          L[, 0:given.d , drop= FALSE]) %*% Q
    }

  
  ## prepare return-values       
  d.seq <- seq.int(0, (max.rk-1)) # dimension-sequence
  E     <- M.Eval[1:max.rk]
  sqr.E <- sqrt(E)
  sum.e <- sum(E)
  cum.e <- cumsum(E)
  V.d   <- c(sum.e, sum.e-cum.e[-length(cum.e)])

  ## return
  structure(list(L              = L,
                 R              = R,
                 Q.orig         = Q.non.smth,
                 Q.orig.smth    = Q,
                 Q.fit          = Q.fit,
                 spar.low       = spar.low,
                 E              = E,
                 sqr.E          = sqr.E,
                 given.d        = given.d,
                 d.seq          = d.seq,
                 V.d            = V.d,
                 nr             = nr,
                 nc             = nc ,
                 cov.mat        = M,
                 dual           = NULL,
                 fpca.method    = "Kneip"),
            class   = "fsvd.pca")
}




############################
# Main function: fpac.fit  #
# Descrition:
# The function calculates functional spectral decomposition of a matrix in >>dat<<
# dimension: dim(dat) == T x N
# and, according to the argument >>given.d<<, does pca-fitting. Generally it is assumend
# that the column of the TxN-matrix >>dat<< are demaend by the "mean-column".
# Calls:                   #
# is.regular.panel()       #
# fsvc.pca()               #
# restrict.pca()           #
# Takes:===========================================================================================#
# dat                      (smoothed raw-data for fPCA)                                            #
# given.d = NULL           (user given dimension)                                                  #
# restrict.mode            ("restrict.factors" or "restrict.loadings")                             #
# allow.dual = TRUE        (possibility to switch of the dual-matrix calculations)                 #
# neglect.neg.ev = TRUE    (if TRUE:  max.rk = number of pos. eigenvalues                          #
#                           if FALSE: max.rk = number of pos. and evtl. numerical neg. eigenvalues)#
# Gives:===========================================================================================#
# factors
# loadings         
# fitted.values    
# orig.values.smth
# orig.values
# spar.low
# cov.matrix      
# eigen.values    
# Sd2              
# given.d        
# data.dim   
# dual       
# L          
####################################################################################################

fpca.fit <- function(dat,
                     fpca.method    = c("Ramsay", "Kneip"),
                     spar.low       = NULL,# no smoothing, iff: spar.low=0  
                     given.d        = NULL,
                     restrict.mode  = c("restrict.factors","restrict.loadings"),
                     allow.dual     = TRUE,
                     neglect.neg.ev = TRUE){

  ## Method: Ramsay or Kneip?
  fpca.method  <- match.arg(fpca.method)  
  ## Check input
  is.regular.panel(dat, stopper = TRUE)

  ## fPCA
  fpca.obj       <- switch(fpca.method,
                           Ramsay  =  fsvd.pca.ramsay(dat,
                                                      given.d        = given.d,
                                                      allow.dual     = allow.dual,
                                                      neglect.neg.ev = neglect.neg.ev),
                           Kneip   = fsvd.pca.kneip(dat,
                                                    given.d          = given.d,
                                                    neglect.neg.ev   = neglect.neg.ev))
  
  ## impose Restrictions (default: restrict.factors such that F'F/T = I )
  result        <- restrict.pca(fpca.obj)

  ## return
  structure(result, class = "fpca.fit")
}
