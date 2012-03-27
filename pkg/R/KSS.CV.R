KSS.CV <- function(kappa.interv, Y, X, N, T, P, tol=tol){
  ## Y (TN x 1)
  ## X (TN x P)
  Outer.CV <- function(kappa){
    Inner.CV <- function(i, kappa){
      Y.mat <- matrix(Y, T,     N)					    # (T x N)
      X.mat <- matrix(X, T, (N*P))                                          # (T x N*P)
      
      Y.mat.min_i <- Y.mat[,-i]                                             # (T x N-1)
      X.mat.min_i <- X.mat[,-seq(from=c(i),to=c(N*P),by=N)]                 # (T x (N-1)*P)
      
      Y.min_i <- matrix(Y.mat.min_i, ncol=1) # (T(N-1) x 1)
      X.min_i <- matrix(X.mat.min_i, ncol=P) # (T(N-1) x P)
      
      Y.mat.min_i.smth  <- smooth.Pspline(x = seq(0,1,len=T), y = Y.mat.min_i, spar = kappa)$ysmth  #(T x (N-1))    
      X.mat.min_i.smth  <- smooth.Pspline(x = seq(0,1,len=T), y = X.mat.min_i, spar = kappa)$ysmth  #(T x (N-1)P)
           
      ## calculate beta coefficents 
      Y.min_i.smth        <- matrix(Y.mat.min_i.smth,  nrow= ((N-1)*T), ncol = 1)      # (T(N-1) x 1)
      X.min_i.smth        <- matrix(X.mat.min_i.smth,  nrow= ((N-1)*T), ncol = P)      # (T(N-1) x P)
           
      t.X.X.min_i         <- crossprod(X.min_i)       			               # (PxP)
      t.X.X.min_i.smth    <- crossprod(X.min_i, X.min_i.smth)		               # (PxP)
           
      t.X.Y.min_i         <- crossprod(X.min_i, Y.min_i)     		               # (Px1)
      t.X.Y.min_i.smth    <- crossprod(X.min_i, Y.min_i.smth)   		       # (Px1)
      
      bloc1               <- t.X.X.min_i - t.X.X.min_i.smth     		       # (PxP)
      bloc2               <- t.X.Y.min_i - t.X.Y.min_i.smth     	               # (Px1)

      ## common-Slope.Coefficients:
      com.slops.0.min_i   <- solve(bloc1)%*%bloc2				           # (Px1)
      ## calculate first step residuals and estimate dimension of factor-structure
      Residu.mat          <- matrix((Y.min_i - (X.min_i %*% com.slops.0.min_i)), T, (N-1)) # (Tx(N-1))        
      ## functional pca
      fpca.fit.obj     <- fpca.fit(Residu.mat, spar=kappa)
      d.hat            <- c(OptDim(Obj=Residu.mat, criteria="KSS.C", spar=1.362e-05)$summary)# 
      print(d.hat)
      Reminder_i       <- Y.mat[,i] - X.mat[,seq(from=c(i),to=c(N*P),by=N)] %*% com.slops.0.min_i      
      if(d.hat == 0){
        Sum.Resid_i    <- sum(Reminder_i^2)
      }else{
        factors        <- fpca.fit.obj$factors[,0:min(d.hat,round(sqrt(min(N, T)))), drop= FALSE]
        Sum.Resid_i    <- sum(residuals(lm(Reminder_i~factors))^2)
      }
    }
    result <-  sum(apply(matrix(1:N, N, 1),1, Inner.CV, kappa=kappa))
    cat(".")
    result
  }
  cat("Progress: CV-Optimization is running.\n")
  return.obj <- optimize(f=Outer.CV, interval=kappa.interv, tol=tol)
  return(return.obj) 
}



## KSS.CV <- function(kappa.interv, Y, X, N, T, P){
  
##   ## Y (TN x 1)
##   ## X (TN x P)
##   Inner.CV <- function(kappa){
##     Y.mat             <- matrix(Y, T,     N)# (T x N)
##     X.mat             <- matrix(X, T, (N*P))# (T x N*P)

##     Y.mat.min_i.array <- array(data= NA, dim=c(T,(N-1),    N))# (T x (N-1)     x N)
##     X.mat.min_i.array <- array(data= NA, dim=c(T,((N-1)*P),N))# (T x ((N-1)*P) x N)

##     Y.min_i.array     <- array(data= NA, dim=c(T*(N-1), 1, N))# (T(N-1) x 1 x N)
##     X.min_i.array     <- array(data= NA, dim=c(T*(N-1), P, N))# (T(N-1) x P x N)
##     for(i in 1:N){
##       Y.mat.min_i.array[,,i] <- Y.mat[,-i]                            # (T x (N-1)   x N)  				    
##       X.mat.min_i.array[,,i] <- X.mat[,-seq(from=c(i),to=c(N*P),by=N)]# (T x (N-1)*P x N)
      
##       Y.min_i.array[,,i]     <- matrix(Y.mat[,-i], (T*(N-1)), 1)                            # (T(N-1) x 1 x N)
##       X.min_i.array[,,i]     <- matrix(X.mat[,-seq(from=c(i),to=c(N*P),by=N)],(T*(N-1)), P) # (T(N-1) x P x N)
##     }
    
##     Y.mat.min_i.smth <- array(NA, dim=c(T, (N-1),     N))
##     X.mat.min_i.smth <- array(NA, dim=c(T, ((N-1)*P), N))
##     for(i in 1:N){
##       Y.mat.min_i.smth[,,i]<-smooth.Pspline(x=seq.int(1,T),y=Y.mat.min_i.array[,,i],spar=kappa)$ysmth #(T x (N-1)   x N) 
##       X.mat.min_i.smth[,,i]<-smooth.Pspline(x=seq.int(1,T),y=X.mat.min_i.array[,,i],spar=kappa)$ysmth #(T x (N-1)*P x N)
##     }
      
##     ## calculate beta coefficents
##     Y.min_i.smth.array <- array(NA, dim=c((T*(N-1)), 1, N))
##     X.min_i.smth.array <- array(NA, dim=c((T*(N-1)), P, N))
##     for(i in 1:N){
##       Y.min_i.smth.array[,,i] <- matrix(Y.mat.min_i.smth[,,i],  nrow= ((N-1)*T), ncol = 1)# (T(N-1) x 1 x N)
##       X.min_i.smth.array[,,i] <- matrix(X.mat.min_i.smth[,,i],  nrow= ((N-1)*T), ncol = P)# (T(N-1) x P x N)
##     }

##     t.X.X.min_i.array      <- array(NA, dim=c(P,P,N))
##     t.X.X.min_i.smth.array <- array(NA, dim=c(P,P,N))
##     t.X.Y.min_i.array      <- array(NA, dim=c(P,1,N))
##     t.X.Y.min_i.smth.array <- array(NA, dim=c(P,1,N))

##     bloc1.array            <- array(NA, dim=c(P,P,N))
##     bloc2.array            <- array(NA, dim=c(P,1,N))
##     for(i in 1:N){
##       t.X.X.min_i.array[,,i]         <- crossprod(X.min_i.array[,,i])       		      # (PxPxN)
##       t.X.X.min_i.smth.array[,,i]    <- crossprod(X.min_i.array[,,i], X.min_i.smth.array[,,i])# (PxPxN)
      
##       t.X.Y.min_i.array[,,i]         <- crossprod(X.min_i.array[,,i], Y.min_i.array[,,i])     # (Px1xN)
##       t.X.Y.min_i.smth.array[,,i]    <- crossprod(X.min_i.array[,,i], Y.min_i.smth.array[,,i])# (Px1xN)
      
##       bloc1.array[,,i]               <- t.X.X.min_i.array[,,i] - t.X.X.min_i.smth.array[,,i]  # (PxPxN)
##       bloc2.array[,,i]               <- t.X.Y.min_i.array[,,i] - t.X.Y.min_i.smth.array[,,i]  # (Px1xN)
##     }
##     ## common-Slope.Coefficients:
##     com.slops.0.min_i.mat <- matrix(NA, P, N)
##     for(i in 1:N){
##       com.slops.0.min_i.mat[,i]   <- solve(bloc1.array[,,i])%*%bloc2.array[,,i]		 # (PxN)
##     }
##     ## calculate first step residuals and estimate dimension of factor-structure
##     Residu.array <- array(NA, dim=c(T, (N-1), N))
##     for(i in 1:N){
##       Residu.array[,,i]<-matrix((Y.min_i.array[,,i]-(X.min_i.array[,,i]%*%com.slops.0.min_i.mat[,i])),T,(N-1))#(Tx(N-1)xN)
##     }        
##       ## functional pca
##     fpca.fit.obj <- vector("list",N)
##     for(i in 1:N){
##       fpca.fit.obj[[i]] <- fpca.fit(Residu.array[,,i], spar=kappa)
##     }
##     d.hat <- NULL
##     for(i in 1:N){
##       d.hat          <- c(d.hat,c(OptDim(Obj=Residu.array[,,i], criteria="KSS.C", spar=0)$summary))#Vector-length N
##     }
##     Reminder_i.mat <- matrix(NA, T,N)
##     for(i in 1:N){
##       Reminder_i.mat[,i]  <- Y.mat[,i] - X.mat[,seq(from=c(i),to=c(N*P),by=N)]%*%com.slops.0.min_i.mat[,i]
##     }
##     Sum.Resid_i <- numeric(N)
##     for(i in 1:N){
##       if(d.hat[i] == 0){
##         Sum.Resid_i[i]    <- sum(Reminder_i.mat[,i]^2)
##       }else{
##         factors           <- fpca.fit.obj[[i]]$factors[,0:d.hat[i], drop= FALSE]
##         Sum.Resid_i[i]    <- sum(residuals(lm(Reminder_i.mat[,i]~factors))^2)
##       }
##     }
##     result <- sum(Sum.Resid_i)
##     result
##   }
##   return.obj <- optimize(f=Inner.CV, interval=kappa.interv)
##   return(return.obj) 
## }



