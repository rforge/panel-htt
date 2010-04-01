pca.fit <-
function(Q, given.d=NULL, calcul.scores = TRUE)
 {
 nr   <- nrow(Q)
 nc   <- ncol(Q)
 dual       <- ifelse(nr>nc, TRUE, FALSE)
 if(dual) Q <- t(Q)

  # Compute spectral decomposion 

 cov.mat <- tcrossprod(Q)
 Spdec   <- eigen(cov.mat)
 Eval   <- Spdec[[1]]
 Evec   <- Spdec[[2]]

  # compare rank and given.d

 max.rk  <- length(Eval[Eval > 0])
 if(is.null(given.d)) given.d <- max.rk
 else 
  {
  if(given.d > max.rk) warning(c("The given dimension 'given.d' is larger than the number of positve eigen values."))
  given.d <- min(given.d, max.rk)
  }

  # compute spectral variance decomposition

 U      <- Evec[,1:max.rk , drop= FALSE]
 sqr.E  <- sqrt(Eval[1:max.rk ])

 if(calcul.scores)
  {
  S <- crossprod(Q, U)[, 1:given.d , drop= FALSE]
  W <- S %*% diag(sqr.E^{-1}) 
  Q.fit  <- tcrossprod(U, S)
  }
 else
  {
  W <- NULL
  Q.fit  <- NULL
  }

  # convert dimension if dual covariance matrix is  used

 if(dual)
  {
  u <- U
  U <- W
  W <- u
  Q.fit <- t(Q.fit)
  }

  # about dimension
 d.seq <- seq.int(0, (given.d-1))
 E <- Eval[1:given.d]
 sum.e <- sum(E)
 cum.e <- cumsum(E)
 V.d   <- c(sum.e, sum.e-cum.e[-length(cum.e)])

 structure(list(U=U, W=W, Q.fit=Q.fit, E=E, sqr.E=sqr.E, given.d=given.d
  , d.seq=d.seq, V.d=V.d, nr=nr, nc=nc ,cov.mat=cov.mat, dual=dual)
  , class = "pca.fit")
      }

