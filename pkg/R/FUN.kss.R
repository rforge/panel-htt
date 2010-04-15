FUN.kss <- function(formula, given.d, x.weights=NULL, norder=2, df=norder+2, spar=0, method=1 
	, effect = c("pooled", "individual", "time", "twoways")
	,... ) # ...  = sign.vec= NULL
	{

  # extract data from formula
	obj.Pformula <- FUN.Pformula(formula= formula, effect = effect)
        
	N <- ncol(obj.Pformula[[1]]$ODM)
	T <- nrow(obj.Pformula[[1]]$ODM)
	P <- length(obj.Pformula)-1
        TRdat <- sapply(1:(P+1), function(i) obj.Pformula[[i]]$TDM) 

	TR.Y  <- TRdat[, 1, drop = FALSE]						# (TN x 1)
	TR.X  <- TRdat[, 2:(P+1), drop = FALSE]					        # (TN x P)

	TR.Y.mat <- matrix(TR.Y, T, N)						        # (T x N)
	TR.X.mat <- matrix(TR.X, T, (N*P))						# (T x NP)
	
  # smooth.splines with undersmoothing
#########################################################################################################
  # undersmoothing crit fehlt noch daher erstmal noch dr=10    
        
        TR.Y.mat.smth  <- smooth.Pspline(x = seq.int(1,T) 
						, y = TR.Y.mat				#(T x N)
						, df=10, method=2)$ysmth
        TR.X.mat.smth  <- smooth.Pspline(x = seq.int(1,T) 
						, y = TR.X.mat
                                         , df=10, method=2)$ysmth                       #(T x NP)
#########################################################################################################
  # calculate beta coefficents

        TR.Y.smth<- matrix(TR.Y.mat.smth,     nrow= (N*T), ncol = 1)	               # (TN x 1)
        TR.X.smth<- matrix(TR.X.mat.smth, nrow= (N*T), ncol = P)	               # (TN x P)

        t.TR.X.TR.X      <- crossprod(TR.X)       			               # (PxP)
        t.TR.X.TR.X.smth <- crossprod(TR.X, TR.X.smth)		        	       # (PxP)

        t.TR.X.TR.Y      <- crossprod(TR.X, TR.Y)     		                       # (Px1)
        t.TR.X.TR.Y.smth <- crossprod(TR.X, TR.Y.smth)   		               # (Px1)

        bloc1 <- t.TR.X.TR.X - t.TR.X.TR.X.smth     			               # (PxP)
        bloc2 <- t.TR.X.TR.Y - t.TR.X.TR.Y.smth     		                       # (Px1)

        com.slops.0 <- solve(bloc1)%*%bloc2					       # (Px1)
        #com.slops.0 <- matrix(com.slops.0, nrow= P,ncol = 1)

  # calculate first step residuals

	Residu.mat    <- matrix((TR.Y - TR.X %*% com.slops.0), T, N)

        obj.gafsa <- gafsa(dat = Residu.mat, dim.crit="KSS.C")



  # additive effects      
  

        
      }
