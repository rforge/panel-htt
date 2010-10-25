########### spectral variance decomposition and pca fitting if asked ###########################

fsvd.pca <- function(Q,  
                     given.d         = NULL,
                     calcul.loadings = TRUE,
                     allow.dual      = TRUE,
                     neglect.neg.ev  = FALSE,
                     spar.low        = NULL # There would be no smoothing, iff: spar.low=0
                     ){

  ## extract data information
  
  nr      <- nrow(Q)
  nc      <- ncol(Q)
  if(nr>nc && calcul.loadings && allow.dual){
    dual  <- TRUE
  }else{
    dual  <- FALSE
  }

  ## save original Q-values 

  Q.non.smth <- Q 

  ## smoothing Q (small degree of undersmoothing)==============================================#
  if(is.null(spar.low)){                                                                       #
    spar.low <- smooth.Pspline(x=seq(0, 1, length.out=nr), y=Q, method = 3       )$spar  * 0.8 #
  }                                                                                            #
  Q          <- smooth.Pspline(x=seq(0, 1, length.out=nr), y=Q, spar   = spar.low)$ysmth       #
  ##===========================================================================================#

  ## do dual-transformation
  
  if(dual){
    Q          <- t(Q)
    # damit die daten bei der rückgabe die gleichen dimensionen haben
    Q.non.smth <- t(Q.non.smth)
         }
    
  ## trapezoidal rule (for integral-approximation)=============================#
  ## hier sehr einfach, da: len.Interval==n.discr-1                            #
  ## Achtung hier wird angenommen, dass beide indizes i bzw. t bei 1 beginnen! #
  ## daher nc-1 bzw. nr-1 in für len.Interval                                  #
  len.Interval    <- ifelse(dual, nc-1, nr-1)                                  #
  n.discr         <- ifelse(dual, nc,   nr  )                                  #  
  h               <- (len.Interval)/(n.discr-1)      # hier: h==1              #
  w               <- c(h/2, rep(h, n.discr-2), h/2)  # hier: w=1/2,1,...,1,1/2 #
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
  L                         <- Evec[,1:max.rk , drop= FALSE]
  ## *fun*ctional approximation
  L.fun                     <- diag(1/(sqrt(w))) %*% L

  ## sqrt of eigenvalues      
  if(!neglect.neg.ev){
    sqr.E                   <- c(sqrt(Eval[Eval > 0]), rep(0, (max.rk - nbr.pos.ev)))
  }else{
    sqr.E	            <- sqrt(Eval[1:max.rk ])
  }

  ## loadings-parameter (scores)
  if(calcul.loadings){
    S                       <- crossprod(Q, L)[, 0:max.rk , drop= FALSE]     # crossprod: t(Q) %*% L
    ## scaling such that: ||R.fun||_E == I
    R                       <- S %*% diag(1/(sqrt(diag(crossprod(S)))))
    ## approximation to L2-norm
    R.fun                   <- diag(1/(sqrt(w))) %*% R

    ## no pca-fitting
    if(((given.d==max.rk)&&!neglect.neg.ev)){
      Q.fit                 <- Q
    ## pca-fitting  
    }else{
      Q.fit                 <- tcrossprod(L.fun[, 0:given.d , drop= FALSE],
                                          L.fun[, 0:given.d , drop= FALSE]) %*% Q
    }
  ## no loadings-parameter (scores)  
  }else{
    R.fun                   <- NULL
    if(((given.d==max.rk)&&!neglect.neg.ev)){
      Q.fit <- Q
    }else{
      Q.fit                 <- tcrossprod(L.fun[, 0:given.d , drop= FALSE],
                                          L.fun[, 0:given.d , drop= FALSE]) %*% Q
    }
  }

  ## re-convert dimension iff dual covariance matrix was used

  if(dual){
    u          <- L.fun
    L.fun      <- R.fun
    R.fun      <- u
    Q.fit      <- t(Q.fit)
    ## damit die daten bei der rückgabe die gleichen dimensionen haben:
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
                 dual           = dual),
            class   = "fsvd.pca")
}




############################
# Main function: fpac.fit  #
# Descrition:
# The function calculates functional spectral decomposition of a matrix in >>dat<<
# and, according to the argument >>given.d<<, does pca-fitting. If  
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
                     spar.low       = NULL,# no smoothing, iff: spar.low=0  
                     given.d        = NULL,
                     restrict.mode  = c("restrict.factors","restrict.loadings"),
                     allow.dual     = TRUE,
                     neglect.neg.ev = TRUE){
  ## Check input
  is.regular.panel(dat, stopper = TRUE)
  ## fPCA
  fsvd.pca.obj  <- fsvd.pca(dat,
                            given.d        = given.d,
                            allow.dual     = allow.dual,
                            neglect.neg.ev = neglect.neg.ev)
  ## impose Restrictions (default: restrict.factors such that F'F/T = I )
  result        <- restrict.pca(fsvd.pca.obj)

  ## return
  structure(result, class = "fpca.fit")
}
