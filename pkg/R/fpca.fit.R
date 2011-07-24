#########################################################################################
## Ramsay-approach 
#########################################################################################

fsvd.pca <- function(Q,
                            allow.dual      = TRUE,
                            given.d         = NULL,
                            calcul.loadings = TRUE,
                            neglect.neg.ev  = FALSE){
  ## extract data information
  nr      <- nrow(Q)
  nc      <- ncol(Q)

  ## save original Q-values 
  Q.non.smth <- Q 

  ## smoothing Q (small degree of undersmoothing) ===============================================#
  spar.low <- smooth.Pspline(x=seq(0, 1, length.out=nr), y=Q, spar=0.01, method = 4       )$spar  * 0.99    #
                                                                                                 #
  Q          <- smooth.Pspline(x=seq(0, 1, length.out=nr), y=Q, spar   = spar.low)$ysmth         #
  ##=============================================================================================#
  
  ## For method=Ramsay: dual-matrix or not?
  dual    <- (nr>nc & allow.dual | calcul.loadings )
 
  if(dual){
    Q <- t(Q)
  }

##   ## classical approach of Ramsay et al.
##   ## trapezoidal rule (for integral-approximation)=============================#
##   ## Hier nicht unbedingt notwendig, da                                        #
##   ## angenommen wird, dass:                                                    #
##   ## -beide indizes i bzw. t bei 1 beginnen!                                   #
##   ## -len.Interval==n.discr-1 sodass h=(len.Interval)/(n.discr-1) = 1          #
##   len.Interval    <- ifelse(dual, nc-1, nr-1)       # dual: N-1, non-dual: T-1 #
##   n.discr         <- ifelse(dual, nc,   nr  )       # dual: N  , non-dual: T   #  
##   h               <- (len.Interval)/(n.discr-1)                                #
##   w               <- c(rep(h, n.discr))   #c(h/2, rep(h, n.discr-2), h/2)              
                                        #
##  cov.mat         <- diag(sqrt(w)) %*% tcrossprod(Q) %*% diag(sqrt(w))         #
    cov.mat         <- tcrossprod(Q)
  ##===========================================================================#
  
  ## Compute spectral decompostion 
  Spdec           <- eigen(cov.mat, symmetric= TRUE)
  Eval	          <- Spdec[[1]]
  Evec    	  <- Spdec[[2]]

  ## 

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
##  L.fun                     <- diag(1/(sqrt(w))) %*% L
  L.fun                     <- L

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
    ## scaling such that: ||R.fun||_E == 1
    R                       <- S %*% diag(1/(sqrt(diag(crossprod(S)))))  # if dual: dim(R)= TxN, if non-dual: dim(R)=NxT
##     ## approximation to L2-norm ===========================================================================#
##     len.Interval.ast    <- ifelse(!dual, nc-1, nr-1)                     # if dual: N-1, if non-dual: T-1  #
##     n.discr.ast         <- ifelse(!dual, nc,   nr  )                     # if dual: N  , if non-dual: T    #
##     h.ast               <- (len.Interval.ast)/(n.discr.ast-1)                                              #
##     w.ast               <- c(rep(h.ast, n.discr.ast)) #c(h.ast/2, rep(h.ast, n.discr.ast-2), h.ast/2)                                  #
##     R.fun               <- diag(1/(sqrt(w.ast))) %*% R                                                     #
    R.fun               <- R                 
##     ##=====================================================================================================#

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
    ## no pca-fitting
    if(((given.d==max.rk)&&!neglect.neg.ev)){
      Q.fit <- Q
    ## pca-fitting
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
                     given.d        = NULL,
                     restrict.mode  = c("restrict.factors","restrict.loadings"),
                     allow.dual     = TRUE,
                     neglect.neg.ev = TRUE){


  ## Check input
  is.regular.panel(dat, stopper = TRUE)

  ## fPCA
  fpca.obj       <- fsvd.pca(Q              = dat,
                             given.d        = given.d,
                             allow.dual     = allow.dual,
                             neglect.neg.ev = neglect.neg.ev)
                           
  
  ## impose Restrictions (default: restrict.factors such that F'F/T = I )
  result        <- restrict.pca(fpca.obj)

  ## return
  structure(result, class = "fpca.fit")
}
