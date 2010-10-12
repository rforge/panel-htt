
################################## Entire Updated Estimation ########################################
# Input:
#	  1. dat.matrix	=  is a matrix where the first colomn containing the (NTx1) vector of Ys
#		  and the remaining colomns containing the (NTxP) vector of Xs
# 	  2. dat.dim         dat.dim[1] is nr (nbr of rows) and dat.dim[2] is nc (nbr of colomns) 
#          of the panel matrix Y.
#	  3. dim.criterion	= c("PC1", "PC2", "PC3", "IC1", "IC2" , "IC3", "IPC1", "IPC2", "IPC3"
#	     , "KSS.C1", "KSS.C2", "ED", "ER", "GR")
#	  4. factor.dim	= the number of factors if it is known (standard is NULL)	
#	  5. d.max 		= the maximum number of factos (an argument needed for the dimension selction)
#	  6. sig2.hat 	= is an argument needed for the dimension selction
#       7. level		= is an argument needed for the dimension selction
#	  8. spar		= is an argument needed for the dimension selction
#	  9. double.iteration =  logical argument, if TRUE the function will estimate the optimal dimension
#	     in an outer iteration the convergence of all the paramters is obtained in a double iteration.
#	     If FALSE the dimension will be estimated parallelly with beta, lambda and F in order 
#	     to reduce the number of computation (desadvantege: convergence to a local a local optimum).
#	     this argument will be neglected if factor.dim is specified.
# Output:
#	  1. $PCA   		   = svd.pca object calculated at the optimal dimension (or given factor.dim)
#	  2. $beta  		   = optimal dimension according to the dimensionalty criterion
#	  3. $opt.d = the slope estimator of the observed regressors (beta.eup)
#	  4. $nbr.iterations = number of iteration
#####################################################################################################

FUN.Eup <- function(dat.matrix, dat.dim, double.iteration = double.iteration, dim.criterion
					, factor.dim, d.max, sig2.hat, level, spar, start.beta) 
  {
#### data 

	y 	<- dat.matrix[, 1, drop = FALSE]
	x 	<- dat.matrix[,-1, drop = FALSE]
	nr 	<- dat.dim[1]
	nc	<- dat.dim[2]
	P	<- ncol(x)

#### start value #####
	beta.0 <- start.beta
	if(is.null(beta.0)) beta.0 <- coef(lm(y~x-1))

#### calculate the inverse once in order to reduce computations of the iterated slope estimator

	FUN.ols.beta <- function(updated.y, x, inv.xx.x){
		beta <- tcrossprod(inv.xx.x, t(updated.y))
		}
	xx 	   <- crossprod(x)
	inv.xx   <- solve(xx)
	inv.xx.x  <- inv.xx%*%t(x)

#### if double iteration factor.dim = d.max

	if(is.null(d.max)) d.max <- round(sqrt(min(nr, nc)))
	given.d <- factor.dim
	if(is.null(factor.dim)) factor.dim <- 1

#### Inner Iteration

	inner.iteration <- function(y, x, inv.xx.x =inv.xx.x , beta.0 = beta.0
				    , factor.dim = factor.dim, d.max=d.max
				    , sig2.hat= sig2.hat, level = level, spar = spar
				    , double.iteration = double.iteration, i=1){
 	#(0)OLS.0: w.0 = y-x%*%beta.0 
		w.0 <- y - tcrossprod(x, t(beta.0))

  	#(0)PCA.0: write w.0 in a matrix form 
		W.0 <- matrix(w.0, nr, nc)

 	#(0)PCA.0: PCA.0 computation, OptDim.0 and y.fitted
		PCA.0    <- svd.pca(W.0, given.d=factor.dim)
		if(!double.iteration && is.null(given.d)){
			OptDim.0   <- EstDim(PCA.0, dim.criterion = dim.criterion 
						, d.max = d.max, sig2.hat= sig2.hat
					 	, level = level, spar = spar)
			opt.dim.0  <- OptDim.0[,2]
			factor.dim <- opt.dim.0
			y.fitted.0 <- PCA.0$L[, 1:opt.dim.0
					, drop = FALSE]%*%diag(PCA.0$sqr.E[1:opt.dim.0]
					, opt.dim.0)%*%t(PCA.0$R[, 1:opt.dim.0, drop = FALSE])	
		}

		else y.fitted.0 = PCA.0$Q.fit

  	#(1)OLS.1: updat y.updated.0 = y - fs.0
		y.updated.1 <- y -  c(y.fitted.0)

  	#(1)OLS.1: OLS.1 computation
		beta.1 <- FUN.ols.beta(y.updated.1, x, inv.xx.x) 

  	# convergence condition
		if(isTRUE(all.equal(beta.0, beta.1, check.attributes = FALSE
			, check.names = FALSE))| i ==150){
			if(i==150) warning(expression("the maximal number of iterations is achieved (150)"))
			Result <- list(PCA=PCA.0, beta=beta.1, factor.dim=factor.dim
			, Nbr.Iterations = i)
			Result
			}
		else inner.iteration(y=y, x=x, inv.xx.x =inv.xx.x , beta.0 = beta.1, 
				    factor.dim = factor.dim, d.max=d.max, 
				    sig2.hat= sig2.hat, level = level, spar = spar,
				    double.iteration = double.iteration, (i+1))
	}

####  outer iterateion 

	entire.iteration <- function(y=y, x=x, inv.xx.x =inv.xx.x , beta.0 = beta.0 
				    , factor.dim = factor.dim, d.max=d.max 
				    , sig2.hat= sig2.hat, level = level, spar = spar
				    , double.iteration = double.iteration
				    , past.iterations = 0, l =1){
	# first inner iteration 
		In.Iter.0 	  <- inner.iteration(y=y, x=x, inv.xx.x =inv.xx.x 
					, beta.0 = beta.0, factor.dim = factor.dim
				      , d.max=d.max, sig2.hat= sig2.hat
				      , level = level, spar = spar
				    	, double.iteration = double.iteration, i =1)
		pca.d0 	  <- In.Iter.0$PCA
		beta.d0 	  <- In.Iter.0$beta 
		opt.d0	  <- In.Iter.0$factor.dim
		nbr.iteration <- In.Iter.0$Nbr.Iterations + past.iterations


		# if double.iteration is TRUE select new opt.d iteratively
		if(double.iteration && is.null(given.d)){	
		  # 1 new optimal dimension
		opt.dim1        <- EstDim(pca.d0, dim.criterion = dim.criterion 
						, d.max = d.max, sig2.hat= sig2.hat
					 	, level = level, spar = spar)
		opt.d1 	    <- opt.dim1[,2]
		  # convergence condition
			if(opt.d1==opt.d0|l == 10){
				if(l==10) warning(expression("the maximal number of outer iterations is achieved (10)"))
				Result  <- list(PCA = pca.d0, beta= beta.d0
						,opt.d =opt.d1, nbr.iterations= nbr.iteration)
				Result
				}
			else entire.iteration(y=y, x=x, inv.xx.x =inv.xx.x , beta.0 = beta.d0
						   , factor.dim = opt.d1, d.max=d.max
						   , sig2.hat= sig2.hat, level = level, spar = spar
						   , double.iteration = double.iteration
						   , past.iterations = nbr.iteration, l = (l+1))
		}
		else {
		Result <- list(PCA = pca.d0, beta= beta.d0
			, opt.d =opt.d0, nbr.iterations=nbr.iteration)
		Result
		}
	}

###### entire iteration result 

	Result	 <- entire.iteration(y=y, x=x, inv.xx.x =inv.xx.x , beta.0 = beta.0 
			    , factor.dim = factor.dim, d.max=d.max 
			    , sig2.hat= sig2.hat, level = level, spar = spar
			    , double.iteration = double.iteration, past.iterations = 0, l =1)

	Result
  }

############################### Eup.default ########################################################
# Input:
#	  1. Formel
#	  2. additive.effects = c("none", "individual", "time", "twoways")
#	  3. dim.criterion	= c("PC1", "PC2", "PC3", "IC1", "IC2" , "IC3", "IPC1", "IPC2", "IPC3"
#	     , "KSS.C1", "KSS.C2", "ED", "ER", "GR")
#	  4. factor.dim	= the number of factors if it is known (standard is NULL)	
#	  5. d.max 		= the maximum number of factos (an argument needed for the dimension selction)
#	  6. sig2.hat 	= is an argument needed for the dimension selction
#       7. level		= is an argument needed for the dimension selction
#	  8. spar		= is an argument needed for the dimension selction
#	  9. double.iteration =  logical argument, if TRUE the function will estimate the optimal dimension
#	     in an outer iteration the convergence of all the paramters is obtained in a double iteration.
#	     If FALSE the dimension will be estimated parallelly with beta, lambda and F in order 
#	     to reduce the number of computation (desadvantege: convergence to a local a local optimum).
#	     this argument will be neglected if factor.dim is specified.
# Output:
#	  1. $PCA   		= svd.pca object calculated at the optimal dimension (or given factor.dim)
#	  2. $beta  		= the slope estimator of the observed regressors (beta.eup)
#	  3. $opt.d			= optimal dimension according to the dimensionalty criterion
#	  4. $nbr.iterations 	= number of iteration
#####################################################################################################

Eup.default <- function(formula, additive.effects = c("none", "individual", "time", "twoways")
				, dim.criterion	= c("PC1", "PC2", "PC3", "IC1", "IC2" , "IC3", "IPC1"
				, "IPC2", "IPC3" , "KSS.C1", "KSS.C2", "ED", "ER", "GR")
				, d.max = NULL,  sig2.hat = NULL, level= 0.01, spar=NULL, factor.dim = NULL
				, double.iteration = FALSE, start.beta= NULL
				, restrict.mode= c("restrict.factors","restrict.loadings")){
### substruct data from fomrmula and perfome a transformation according to additive.effects 
  # check fomula
	if(!class(formula)=="formula"){
		stop("\n Argument >>formula<< needs a formula-object like y~x1+... where the elements are matrices")
	}	
  # set additive.effects	
	additive.effects <- match.arg(additive.effects)
  # transform the data (output in a list)
	PF.obj <- FUN.Pformula(formula = formula, effect = additive.effects)    
	nc 		<- ncol(PF.obj[[1]]$ODM)
	nr 		<- nrow(PF.obj[[1]]$ODM)
	P  		<- length(PF.obj)-1
  # prepare for FUN.default (output in a matrix)
	dat.dim 	  <- c(nr, nc)
	dat.matrix	  <- sapply(1:(P+1), function(i) PF.obj[[i]]$TDM) 
	dim.criterion <- match.arg(dim.criterion)
	
  # Estimation results

	tr.model.est <- FUN.Eup(dat.matrix		= dat.matrix
					, dat.dim		= dat.dim
					, double.iteration= double.iteration
					, start.beta	= start.beta
					, dim.criterion 	= dim.criterion
					, factor.dim	= factor.dim
					, d.max		= d.max
					, sig2.hat		= sig2.hat
					, level		= level
					, spar			= spar) 

    # Eup beta and Nbr.iteration
	Nbr.iteration	<- tr.model.est$nbr.iterations
	beta.Eup		<- tr.model.est$beta

    # factor dimension 
	proposed.dim <- ifelse(is.null(factor.dim), "not indicated", factor.dim)			
	if(is.null(factor.dim)){
		optimal.dim <- tr.model.est$opt.d
		}
	else{
		optimal.dim <- EstDim(tr.model.est$PCA, dim.criterion = dim.criterion 
						, d.max = d.max, sig2.hat= sig2.hat
					 	, level = level, spar = spar)[,2]
#		if(optimal.dim < factor.dim) optimal.dim <- c(" The used dimensionalty criterion "
#						  , "(" , dim.criterion , ")" ," detected a smaller number of factors "
#						  , "(", optimal.dim, ")", "as proposed " , "(" , dim.criterion , ")." 
#						  , " The Eup estimators may be efficienter.")
#		else{
#			if(optimal.dim > factor.dim) optimal.dim  <- c(" The used dimensionalty criterion "
#						  , "(" , dim.criterion , ")" ," detected a larger number of factors "
#						  , "(", optimal.dim, ").", " The estimated parameters may be baised.")
#			else optimal.dim  <- c(" The used dimensionalty criterion "
#						  , "(" , dim.criterion , ")" ," detected a number of factors equal 
#						  to the proposed dimension", "(" , dim.criterion , ").")
#		 	}
		}
	used.dim <- tr.model.est$opt.d


    # factor structure and resuduals 
	restrict.fs.a.resid	<- restrict.pca(tr.model.est$PCA, restrict.mode=restrict.mode)					
	fs.and.resid		<- restrict.fs.a.resid$orig.values
	factors			<- restrict.fs.a.resid$factors
	loadings			<- restrict.fs.a.resid$loadings
	unobs.fact.struct		<- restrict.fs.a.resid$fitted.values
	residuals			<- fs.and.resid - unobs.fact.struct

## Inference about beta

	slope.inf 		<-Eup.inference(dat.matrix, dat.dim, used.dim, beta.Eup, factors, loadings, residuals)

## Results
				
	list(slope.inf, used.dim, proposed.dim, optimal.dim, Nbr.iteration)
  }

################################## Eup slope inference ####################
## Input:
# 		dat.matrix	= the data in matrix form where the first colomn contain NT vector of Y second one the NT vector of the first regressor X
# 		dat.dim	= the dimension of the data N and T
# 		used.d	= used dimension d in the interative procedure
#	 	beta.Eup	= the estimated slope estimator for given d
# 		factors	= common factors after scaling according the used restriction
# 		loadings	= individual loading parameters after scaling (according restriction)
#		residuals	= the residual terms
## Output:
#		Eup slope Estimate
#		std
#		Pr(>|z|)
###########################################################################
Eup.inference <- function(dat.matrix, dat.dim, used.dim, beta.Eup
				, factors, loadings, residuals){
	y 	<- dat.matrix[, 1, drop = FALSE]
	x 	<- dat.matrix[,-1, drop = FALSE]
	nr 	<- dat.dim[1]
	nc	<- dat.dim[2]
	P	<- ncol(x)
	d	<- used.dim
	slope <- beta.Eup
	F	<- factors
	A	<- loadings
	res	<- residuals


## Projection matrix of the factors F
	if(all(d == 0)) P.F = 0
	else{
 		FtFm1 <- diag(crossprod(F))^{-1}
		P.F   <- F%*%diag(FtFm1, d)%*%t(F)
		}
	M.F   <- diag(1, ncol= nr, nrow= nr) - P.F

## Projection matrix of the loadings A
	if(all(d == 0)) P.A = 0
	else{
 		AtAm1 <- diag(crossprod(A))^{-1}
		P.A   <- A%*%diag(AtAm1, length(AtAm1))%*%t(A)
		}
	M.A   <- diag(1, ncol= nc, nrow= nc) - P.A

## write the x matrices in a list: each regressor is written in a list component
	X.mat.list <- NULL
	for(p in 1:P) X.mat.list[[p]] <- matrix(x[,p], nr, nc)

  # Z_i = M * X_i - sum{M * X_k*a_ik}/n
	Z.list	<- sapply (X.mat.list, function(X) M.F %*% X %*% M.A , simplify = FALSE)

  # construct the matrix D= sum Z_i'Z_i/NT
	Z		<- sapply (Z.list, function(Z) c(Z), simplify = TRUE)
	ZZ     <- crossprod(Z)/(nr*nc)
	inv.ZZ <- solve(ZZ)
	sig2.hat <- sum(diag( crossprod(res)))/(nr*nc - (nr+nc)*d - p + 1)

	asy.var <- (inv.ZZ * sig2.hat)/(nr*nc)
	mpp <- sqrt(diag(asy.var))
	test <- slope/mpp 
	p.value <- (1 - pnorm(abs(test)))*2
	result <- cbind(slope, mpp,p.value)
	colnames(result) <- c("Eup slope Estimate", "std", " Pr(>|z|)")
	result	
}










