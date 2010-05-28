########### spectral variance decomposition and pca fitting if asked ###########################

fsvd.pca <- function(Q,
                     given.d         = NULL,
                     calcul.loadings = TRUE,
                     allow.dual      = TRUE,
                     neglect.neg.ev  = FALSE){
  
  nr   <- nrow(Q)
  nc   <- ncol(Q)
  if(nr>nc && calcul.loadings && allow.dual) dual = TRUE
  else dual = FALSE
  if(dual) Q      <- t(Q)

  # trapezoidal rule (for integral-approximation)
  
  # Achtung hier wird angenommen, dass beide indizes i bzw. t bei 1 beginnen!== 
  len.Interval    <- ifelse(dual, nc-1, nr-1)                                 #
  n.discr         <- ifelse(dual, nc,   nr  )                                 #  
  h               <- (len.Interval)/(n.discr-1)                               #
  w               <- c(h/2, rep(h, n.discr-2), h/2)                           #
  cov.mat         <- diag(sqrt(w)) %*% tcrossprod(Q) %*% diag(sqrt(w))        #
  #============================================================================

  
  # Compute spectral decomposion 

  Spdec           <- eigen(cov.mat, symmetric= TRUE)
  Eval	          <- Spdec[[1]]
#if(any(Eval<0)) warning(expression("There are negative eigen values!"))
  Evec    	  <- Spdec[[2]]

  # compare rank and given.d
        
  neglect.neg.ev  <- neglect.neg.ev |(dual&&neglect.neg.ev)
  nbr.pos.ev      <- length(Eval[Eval > 0])
  max.rk          <- ifelse(neglect.neg.ev, nbr.pos.ev, length(Eval))
  if(is.null(given.d)) given.d <- max.rk
  else {
    if(given.d > max.rk){
      warning(c("The given dimension 'given.d' is larger than the number of positve eigen values."))
    }
    given.d <- min(given.d, max.rk)
  }

  # compute spectral variance decomposition

  # Approximation of the L2-normalization-constraint:
  #    ||e.fun_i||_L2 = 1 by
  #  w*||e.fun_i||_E  = 1
  
  L                         <- Evec[,1:max.rk , drop= FALSE]
  L.fun                     <- diag(w^-{0.5}) %*% L

        
  if(!neglect.neg.ev) sqr.E <- c(sqrt(Eval[Eval > 0]),rep(0, (max.rk - nbr.pos.ev)))
  else sqr.E	            <- sqrt(Eval[1:max.rk ])
  
  if(calcul.loadings){
    S.fun                   <- crossprod(Q, L.fun)[, 1:max.rk , drop= FALSE]
    R.fun                   <- S.fun %*% diag(diag(crossprod(S.fun))^-{0.5})
    if(((given.d==max.rk)&&!neglect.neg.ev)) Q.fit <- Q
    else Q.fit              <- tcrossprod(L.fun[, 1:given.d , drop= FALSE],
                                          L.fun[, 1:given.d , drop= FALSE]) %*% Q
  }
  
  else{
    R.fun                                          <- NULL
    if(((given.d==max.rk)&&!neglect.neg.ev)) Q.fit <- Q
    else Q.fit                                     <- tcrossprod(L.fun[, 1:given.d , drop= FALSE],
                                                                 L.fun[, 1:given.d , drop= FALSE]) %*% Q
  }

  # convert dimension if dual covariance matrix is used

  if(dual){
    u     <- L.fun
    L.fun <- R.fun
    R.fun <- u
    Q.fit <- t(Q.fit)
  }
    
  # about dimension
        
  d.seq <- seq.int(0, (max.rk-1))
  E     <- Eval[1:max.rk]
  sum.e <- sum(E)
  cum.e <- cumsum(E)
  V.d   <- c(sum.e, sum.e-cum.e[-length(cum.e)])

  structure(list(L.fun   = L.fun,
                 R.fun   = R.fun,
                 Q.fit   = Q.fit,
                 E       = E,
                 sqr.E   = sqr.E,
                 given.d = given.d,
                 d.seq   = d.seq,
                 V.d     = V.d,
                 nr      = nr,
                 nc      = nc ,
                 cov.mat = cov.mat,
                 dual    = dual),
            class   = "fsvd.pca")
}


##########################
# restrict mode "restrict.factors": F'F/T = I
# restrict mode "restrict.loadings" : Lamd'Lamd/N = I 

frestrict.pca <- function(fsvd.pca.obj,
                         restrict.mode = c("restrict.factors","restrict.loadings")){
  
  if(class(fsvd.pca.obj)!="fsvd.pca") stop(c("The fsvd.pca.obj is not a 'fsvd.pca' object"))
  if(is.null(fsvd.pca.obj$R.fun))     stop(c("Loadings-parameter are missing."))
  
 # fsvd.pca object  
  
  cov.mat           <- fsvd.pca.obj$cov.mat
  fitted.values     <- fsvd.pca.obj$Q.fit
  given.d	    <- fsvd.pca.obj$given.d
  L.fun             <- fsvd.pca.obj$L.fun[, 1:given.d , drop= FALSE]
  R.fun             <- fsvd.pca.obj$R.fun[, 1:given.d , drop= FALSE]
  sqr.E             <- fsvd.pca.obj$sqr.E
  E	            <- fsvd.pca.obj$E
  dual              <- fsvd.pca.obj$dual
  nr                <- nrow(fitted.values)
  nc                <- ncol(fitted.values)
  cov.matrix        <- cov.mat/(nr*nc)
  Sd2               <- fsvd.pca.obj$V.d/(nr*nc)
  Eval	            <- E/(nr*nc)

 # restric factors and loadings
  re.mo <-match.arg(restrict.mode)
  switch(re.mo,
         restrict.factors  ={
           factors.fun   <- L.fun  *  sqrt(nr)
           loadings.fun  <- R.fun %*% diag(sqr.E[1:given.d], given.d)/sqrt(nr)
         },
         restrict.loadings ={
           factors.fun   <- L.fun %*% diag(sqr.E[1:given.d], given.d)/sqrt(nc)
           loadings.fun  <- R.fun  *  sqrt(nc)
         }
         )

  list(factors.fun    = factors.fun,
       loadings.fun   = loadings.fun,
       fitted.values  = fitted.values,
       cov.matrix     = cov.matrix,
       eigen.values   = Eval,
       Sd2            = Sd2,
       given.d        = given.d,
       data.dim       = c(nr, nc),
       dual           = dual)		
}

############################
# Main function: fpac.fit  #
# Calls:                   #
# is.regular.panel()       #
# fsvc.pca()               #
# frestrict.pca()          #
# Takes:                   #
# dat                    (smoothed raw-data for fPCA)
# given.d = NULL         (user given dimension)
# restrict.mode          ("restrict.factors" or "restrict.loadings")
# allow.dual = TRUE      (possibility to switch of the dual-matrix calculations)
# neglect.neg.ev = TRUE  (if TRUE:  max.rk = number of pos. eigenvalues
#                         if FALSE: max.rk = number of pos. and evtl. numerical neg. eigenvalues)
############################

fpca.fit <- function(dat,
                     given.d        = NULL,
                     restrict.mode  = c("restrict.factors","restrict.loadings"),
                     allow.dual     = TRUE,
                     neglect.neg.ev = TRUE){
  # Check input
  is.regular.panel(dat, stopper = TRUE)
  # fPCA
  fsvd.pca.obj  <- fsvd.pca(dat,
                            given.d        = given.d,
                            allow.dual     = allow.dual,
                            neglect.neg.ev = neglect.neg.ev)
  # Impose Restrictions
  result        <- frestrict.pca(fsvd.pca.obj)
  
  structure(result, class = "fpca.fit")
}
