######## check the input data if is a regular panel (matrix or dataframe without NA NaN)

is.regular.panel <- function(dat, comment = FALSE, stop = FALSE){
	if(!is.matrix(dat)&!is.data.frame(dat)) {
		if(stop) stop(expression("The data should be in a matrix forme or has a data.frame class"),  call. = FALSE)
		if(comment) message(expression("The data should be in a matrix forme or has a data.frame class"))
		FALSE
		}
	else{
		if(!is.numeric(dat)){
			if(stop) stop(expression("The data is not numeric"),  call. = FALSE)
			if(comment) message(expression("The data is not numeric"))
			FALSE
			}
		else{
			if(any(is.na(dat))| any(is.nan(dat))){
				if(stop) stop(expression("The data has NA or NaN values"),  call. = FALSE)
				if(comment) message(expression("The data has NA or NaN values"))
				FALSE
				}
			else TRUE
			}
		}

	}
############# test
#test <- matrix(c(1,2,0, 3), 2,2 )
#data.frame(test)
#test[1,1] <- "a" 
is.regular.panel(test)

########### spectral variance decomposition and pca fitting if asked ###########################

svd.pca <- function(Q, given.d=NULL, calcul.loadings = TRUE, allow.dual = TRUE)
	{
	nr   <- nrow(Q)
	nc   <- ncol(Q)
	if(nr>nc && calcul.loadings && allow.dual) dual = TRUE
	else dual = FALSE
	if(dual) Q <- t(Q)

  # Compute spectral decomposion 

	cov.mat <- tcrossprod(Q)
	Spdec   <- eigen(cov.mat)
	Eval	  <- Spdec[[1]]
	#if(any(Eval<0)) warning(expression("There are negative eigen values!"))
	Evec	  <- Spdec[[2]]

  # compare rank and given.d

	max.rk  <- ifelse(dual, length(Eval[Eval > 0]), length(Eval[Eval >= 0]))
	if(is.null(given.d)) given.d <- max.rk
	else 
		{
		if(given.d > max.rk) warning(c("The given dimension 'given.d' is larger than the number of positve eigen values."))
		given.d <- min(given.d, max.rk)
		}

  # compute spectral variance decomposition

	U      <- Evec[,1:max.rk , drop= FALSE]
	sqr.E	 <- sqrt(Eval[1:max.rk ])

	if(calcul.loadings)
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
		, class = "svd.pca")
	}


##########################
# restrict mode "restrict.factors": F'F/T = I
# restrict mode "restrict.loadings" : Lamd'Lamd/N = I 

restrict.pca <- function(svd.pca.obj
		, restrict.mode= c("restrict.factors","restrict.loadings"))
	{
	if(class(svd.pca.obj)!="svd.pca") stopt(c("The svd.pca.obj is not a 'svd.pca' object"))
 # svd.pca object  

	cov.mat       <- svd.pca.obj$cov.mat
	fitted.values <- svd.pca.obj$Q.fit
	nr            <- nrow(fitted.values)
	nc            <- ncol(fitted.values)
	dual          <- svd.pca.obj$dual
	U 	       <- svd.pca.obj$U
	W 	       <- svd.pca.obj$W
	cov.matrix    <- cov.mat/(nr*nc)
	sqr.E         <- svd.pca.obj$sqr.E
	E	       <- svd.pca.obj$E
	Sd2           <- svd.pca.obj$V.d/(nr*nc)
	Eval	       <- E/(nr*nc)

 # restric factors and loadings
	re.mo <-match.arg(restrict.mode)
	switch(re.mo,
			restrict.factors={
				factors <- U*sqrt(nr)
				loadings  <- W%*%diag(sqr.E)/sqrt(nr)
			},
			restrict.loadings ={
				factors <- U%*%diag(sqr.E)/sqrt(nc)
				loadings  <- W*sqrt(nc)
			}
		)

	list(factors= factors, loadings= loadings, fitted.values = fitted.values
	, cov.matrix = cov.matrix, eigen.values= Eval ,Sd2 = Sd2
	, data.dim = c(nr, nc), dual=dual)		
	}



pca.fit <- function(dat, given.d=NULL 
		, restrict.mode= c("restrict.factors","restrict.loadings")
		, allow.dual = TRUE){
	is.regular.panel(dat, stop = TRUE)
	svd.pca.obj  <- svd.pca(dat, given.d = given.d, allow.dual= allow.dual)
	result <- restrict.pca(svd.pca.obj)
	structure(result, class = "pca.fit")
	}



obj <- svd.pca(dat)

ob <- pca.fit(test)
is.matrix(ob)