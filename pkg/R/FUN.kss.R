FUN.kss <- function(formula, given.d,
                    effect = c("pooled", "individual", "time", "twoways"),
                    ...) # ...  = sign.vec= NULL
  {

  # extract data from formula
	obj.Pformula <- FUN.Pformula(formula = formula, effect = effect)
        
	N <- ncol(obj.Pformula[[1]]$ODM)
	T <- nrow(obj.Pformula[[1]]$ODM)
	P <- length(obj.Pformula)-1
        
        TRdat <- sapply(1:(P+1), function(i) obj.Pformula[[i]]$TDM) 

	TR.Y  <- TRdat[, 1, drop = FALSE]						# (TN x 1)
	TR.X  <- TRdat[, 2:(P+1), drop = FALSE]					        # (TN x P)

	TR.Y.mat <- matrix(TR.Y, T, N)						        # (T x N)
	TR.X.mat <- matrix(TR.X, T, (N*P))						# (T x NP)
	
  # smooth.splines with undersmoothing
  # undersmoothing: 0.8*GCV-value

        spar.low       <- smooth.Pspline(x = seq.int(1,T), y = TR.Y.mat, method = 3       )$spar * 0.8
        TR.Y.mat.smth  <- smooth.Pspline(x = seq.int(1,T), y = TR.Y.mat, spar   = spar.low)$ysmth       #(T x N)

        spar.low       <- smooth.Pspline(x = seq.int(1,T), y = TR.X.mat, method = 3       )$spar * 0.8
        TR.X.mat.smth  <- smooth.Pspline(x = seq.int(1,T), y = TR.X.mat, spar   = spar.low)$ysmth       #(T x NP)
        
  # calculate beta coefficents

        TR.Y.smth        <- matrix(TR.Y.mat.smth,     nrow= (N*T), ncol = 1)	       # (TN x 1)
        TR.X.smth        <- matrix(TR.X.mat.smth, nrow= (N*T), ncol = P)	       # (TN x P)

        t.TR.X.TR.X      <- crossprod(TR.X)       			               # (PxP)
        t.TR.X.TR.X.smth <- crossprod(TR.X, TR.X.smth)		        	       # (PxP)

        t.TR.X.TR.Y      <- crossprod(TR.X, TR.Y)     		                       # (Px1)
        t.TR.X.TR.Y.smth <- crossprod(TR.X, TR.Y.smth)   		               # (Px1)

        bloc1            <- t.TR.X.TR.X - t.TR.X.TR.X.smth     		               # (PxP)
        bloc2            <- t.TR.X.TR.Y - t.TR.X.TR.Y.smth     	                       # (Px1)

        com.slops.0 <- solve(bloc1)%*%bloc2					       # (Px1)
       #com.slops.0 <- matrix(com.slops.0, nrow= P,ncol = 1)

  # calculate first step residuals

	Residu.mat    <- matrix((TR.Y - TR.X %*% com.slops.0), T, N)

        fAFactMod.obj <- fAFactMod(dat = Residu.mat, dim.crit="KSS")
        return(fAFactMod.obj)

 # reestimate beta
#       factor.stract <- tcrossprod(obj.gafsa$factors, obj.gafsa$scores)  
# 	NEW.TR.Y.mat <- TR.Y.mat - factor.stract
## 	new.y.vec <- as.vector(new.y.mat)
## 	beta <- qr.solve(X.s, new.y.vec)

## # functional intercept
##         if(y.in.list$T=="twoway"){
##           if(is.intercept){
##             f.intercept.indv      <- y.in.list$OVm + y.mean.twoway.indv.vec - x.all.mean.twoway.indv %*% beta
##             f.intercept.indv.smth <- smooth.Pspline(x = basis.t , y = f.intercept.indv, df=df.2, method=2)$ysmth
##             f.intercept.time      <- y.in.list$OVm + y.mean.twoway.time.vec - x.all.mean.twoway.time %*% beta
##             f.intercept.time.smth <- smooth.Pspline(x = basis.t , y = f.intercept.time, df=df.2, method=2)$ysmth
##           }else{
##             f.intercept.indv      <- y.mean.twoway.indv.vec - x.all.mean.twoway.indv %*% beta
##             f.intercept.indv.smth <- smooth.Pspline(x = basis.t , y = f.intercept.indv, df=df.2, method=2)$ysmth
##             f.intercept.time      <- y.mean.twoway.time.vec - x.all.mean.twoway.time %*% beta
##             f.intercept.time.smth <- smooth.Pspline(x = basis.t , y = f.intercept.time, df=df.2, method=2)$ysmth
##           }          
##         }else{
##           if(is.intercept){
##             f.intercept         <- y.in.list$OVm + y.mean.vec - x.all.mean %*% beta
##             f.intercept.smth    <- smooth.Pspline(x = basis.t , y = f.intercept, df=df.2, method=2)$ysmth
##           }else{
##             f.intercept         <- y.mean.vec - x.all.mean %*% beta
##             f.intercept.smth    <- smooth.Pspline(x = basis.t , y = f.intercept, df=df.2, method=2)$ysmth
##           }
##         }
        
## # fitted values

## 	x.beta   <- matrix(x.all.matrix, N*T, P)%*% beta

##         if(y.in.list$T=="twoway"){
##           if(is.intercept){
##             y.fitted <- matrix(y.in.list$OVm, T, N) + matrix(f.intercept.time.smth, T, N) +
##               matrix(f.intercept.indiv.smth, T, N) +  matrix(x.beta, T, N) + factor.stract
##           }else{
##             y.fitted <- matrix(f.intercept.time.smth, T, N) + matrix(f.intercept.indiv.smth, T, N) +
##               matrix(x.beta, T, N) + factor.stract
##           }
##         }else{
##           if(is.intercept){
##             y.fitted <- matrix(y.in.list$OVm, T, N) + matrix(f.intercept.smth, T, N) +
##               matrix(x.beta, T, N) + factor.stract
##           }else{
##             y.fitted <- matrix(f.intercept.smth, T, N) + matrix(x.beta, T, N) + factor.stract
##           }
##         }
                  
##         obj <- list(basis.t, y.matrix, x.all.matrix, N,T,P, df.1, df.2, dim.bf, lambda, eig.func, Residu.mat, Resi.smth.mat)
        
        
## par(mfrow= c(2,3))
## matplot(scores, typ = "l")
## matplot(eig.func, typ = "l")
## matplot(factor.stract, typ= "l")
## plot(f.intercept.smth, typ = "l")
## matplot(y.fitted, typ = "l")
## matplot(y.matrix, typ= "l")
## print(beta)# from simulation: beta = c(-1.5, 5, -1.8)
##         obj


        
  # additive effects      
  

        
      }
