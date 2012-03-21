KSS.CV <- function(kappa.interv, Y, X, N, T, P){
 
  #Y (TN x 1)
  #X (TN x P)
  Outer.CV <- function(kappa){
    Inner.CV <- function(i, kappa){
      Y.mat <- matrix(Y, T,     N)					    # (T x N)
      X.mat <- matrix(X, T, (N*P))                                          # (T x N*P)
      
      Y.mat.min_i <- Y.mat[,-i]                                             # (T x N-1)
      X.mat.min_i <- X.mat[,-seq(from=c(i),to=c(N*P),by=N)]                 # (T x (N-1)*P)
      
      Y.min_i <- matrix(Y.mat.min_i, ncol=1) # (T(N-1) x 1)
      X.min_i <- matrix(X.mat.min_i, ncol=P) # (T(N-1) x P)
      
      Y.mat.min_i.smth  <- smooth.Pspline(x = seq.int(1,T), y = Y.mat.min_i,         spar   = kappa)$ysmth  #(T x (N-1))    
      X.mat.min_i.smth  <- smooth.Pspline(x = seq.int(1,T), y = X.mat.min_i,         spar   = kappa)$ysmth  #(T x (N-1)P)
      X.mat.min_i.smth2 <- smooth.Pspline(x = seq.int(1,T), y = X.mat.min_i.smth,    spar   = kappa)$ysmth  #(T x (N-1)P)
      
      ## calculate beta coefficents 
      Y.min_i.smth        <- matrix(Y.mat.min_i.smth,  nrow= ((N-1)*T), ncol = 1)      # (T(N-1) x 1)
      X.min_i.smth        <- matrix(X.mat.min_i.smth,  nrow= ((N-1)*T), ncol = P)      # (T(N-1) x P)
      X.min_i.smth2       <- matrix(X.mat.min_i.smth2, nrow= ((N-1)*T), ncol = P)      # (T(N-1) x P)
      
      t.X.X.min_i         <- crossprod(X.min_i)       			               # (PxP)
      t.X.X.min_i.smth    <- crossprod(X.min_i, X.min_i.smth)		               # (PxP)
      t.X.X.min_i.smth2   <- crossprod(X.min_i, X.min_i.smth2)		               # (PxP)
      
      t.X.Y.min_i         <- crossprod(X.min_i, Y.min_i)     		               # (Px1)
      t.X.Y.min_i.smth    <- crossprod(X.min_i, Y.min_i.smth)   		       # (Px1)
    
      bloc1               <- t.X.X.min_i - t.X.X.min_i.smth     		       # (PxP)
      bloc2               <- t.X.Y.min_i - t.X.Y.min_i.smth     	               # (Px1)
      ## common-Slope.Coefficients:
      com.slops.0.min_i   <- solve(bloc1)%*%bloc2				       # (Px1)
      ## calculate first step residuals and estimate dimension of factor-structure
      Residu.mat          <- matrix((Y.min_i - (X.min_i %*% com.slops.0)), T, (N-1))   # (Tx(N-1))
    
      ## functional pca
      fpca.fit.obj     <- fpca.fit(Residu.mat, spar=kappa)
      d.hat            <- c(OptDim(Obj=Residu.mat, criteria="KSS.C", spar=kappa)$summary)
      factors          <- fpca.fit.obj$factors[,  1:d.hat, drop= FALSE]
      Reminder_i       <- Y.mat[,i] - X.mat[,seq(from=c(i),to=c(N*P),by=N)] %*% com.slops.0.min_i
      Sum.Resid_i      <- sum(residuals(lm(Reminder_i~factors)))
    }

    sum(sapply(1:N, fun=Inner.CV, kappa=kappa))
  }
  return.obj <- optimize(f=Outer.CV, interval=kappa.interv)
  return(return.obj)
}
