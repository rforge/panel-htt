######## check the input data if is a regular panel (matrix or dataframe without NA NaN)
is.regular.panel <- function(dat,
                             comment = FALSE,
                             stopper = FALSE){
  
  if(!is.matrix(dat)&!is.data.frame(dat)){
    if(stopper) stop(expression("The data should be in a matrix forme or has a data.frame class"),  call. = FALSE)
    if(comment) message(expression("The data should be in a matrix forme or has a data.frame class"))
    FALSE
  }
  else{
    if(!is.numeric(dat)){
      if(stopper) stop(expression("The data is not numeric"),  call. = FALSE)
      if(comment) message(expression("The data is not numeric"))
      FALSE
    }
    else{
      if(any(is.na(dat))| any(is.nan(dat))){
        if(stopper) stop(expression("The data has NA or NaN values"),  call. = FALSE)
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
#is.regular.panel(test)

########### spectral variance decomposition and pca fitting if asked ###########################

svd.pca <- function(Q, given.d=NULL, calcul.loadings = TRUE, allow.dual = TRUE
		, neglect.neg.ev = FALSE)
	{
	nr   <- nrow(Q)
	nc   <- ncol(Q)
	if(nr>nc && calcul.loadings && allow.dual) dual = TRUE
	else dual = FALSE
	if(dual) Q <- t(Q)

  # Compute spectral decomposion 

	cov.mat <- tcrossprod(Q)
	Spdec   <- eigen(cov.mat, symmetric= TRUE)
	Eval	  <- Spdec[[1]]
	#if(any(Eval<0)) warning(expression("There are negative eigen values!"))
	Evec	  <- Spdec[[2]]

  # compare rank and given.d
	neglect.neg.ev <- neglect.neg.ev |(dual&&neglect.neg.ev)
	nbr.pos.ev <- length(Eval[Eval > 0])
	max.rk  <- ifelse(neglect.neg.ev, nbr.pos.ev, length(Eval))
	if(is.null(given.d)) given.d <- max.rk
	else {
		if(given.d > max.rk) warning(c("The given dimension 'given.d' is larger than the number of positve eigen values."))
		given.d <- min(given.d, max.rk)
		}

  # compute spectral variance decomposition

	L      <- Evec[,1:max.rk , drop= FALSE]
	if(!neglect.neg.ev) sqr.E <- c(sqrt(Eval[Eval > 0]),rep(0, (max.rk - nbr.pos.ev)))
	else sqr.E	 <- sqrt(Eval[1:max.rk ])

	if(calcul.loadings){
		S <- crossprod(Q, L)[, 1:max.rk , drop= FALSE]
		R <- S %*% diag(diag(crossprod(S))^-{0.5})
		if(((given.d==max.rk)&&!neglect.neg.ev)) Q.fit <- Q
		else Q.fit  <- tcrossprod(L[, 1:given.d , drop= FALSE]
				, L[, 1:given.d , drop= FALSE])%*%Q
		}

	else{
		R <- NULL
		if(((given.d==max.rk)&&!neglect.neg.ev)) Q.fit <- Q
		else Q.fit  <- tcrossprod(L[, 1:given.d , drop= FALSE]
				, L[, 1:given.d , drop= FALSE])%*%Q
		}

  # convert dimension if dual covariance matrix is  used

	if(dual){
		u <- L
		L <- R
		R <- u
		Q.fit <- t(Q.fit)
		}

  # about dimension
	d.seq <- seq.int(0, (max.rk-1))
	E 	<- Eval[1:max.rk]
	sum.e <- sum(E)
	cum.e <- cumsum(E)
	V.d   <- c(sum.e, sum.e-cum.e[-length(cum.e)])

	structure(list(L=L, R=R, Q.fit=Q.fit, E=E, sqr.E=sqr.E, given.d=given.d
		, d.seq=d.seq, V.d=V.d, nr=nr, nc=nc ,cov.mat=cov.mat, dual=dual)
		, class = "svd.pca")
	}


##########################
# restrict mode "restrict.factors": F'F/T = I
# restrict mode "restrict.loadings" : Lamd'Lamd/N = I 

restrict.pca <- function(svd.pca.obj
		, restrict.mode= c("restrict.factors","restrict.loadings")){
  if(class(svd.pca.obj)!="svd.pca"&&class(svd.pca.obj)!="fsvd.pca"){
    stop(c("The svd.pca.obj is not a 'svd.pca' or a 'svd.pca' object"))
  }
  if(is.null(svd.pca.obj$R))  stop(c("Loadings-parameter are missing."))
        
 # svd.pca object  

	cov.mat       <- svd.pca.obj$cov.mat
	fitted.values <- svd.pca.obj$Q.fit
	given.d	      <- svd.pca.obj$given.d
	L 	        <- svd.pca.obj$L[, 1:given.d , drop= FALSE]
	R 	        <- svd.pca.obj$R[, 1:given.d , drop= FALSE]
	sqr.E         <- svd.pca.obj$sqr.E
	E	        <- svd.pca.obj$E
	dual          <- svd.pca.obj$dual
	nr            <- nrow(fitted.values)
	nc            <- ncol(fitted.values)
	cov.matrix    <- cov.mat/(nr*nc)
	Sd2           <- svd.pca.obj$V.d/(nr*nc)
	Eval	        <- E/(nr*nc)

 # restric factors and loadings
	re.mo <-match.arg(restrict.mode)
	switch(re.mo,
			restrict.factors={
				factors <- L*sqrt(nr)
				loadings  <- R%*%diag(sqr.E[1:given.d], given.d)/sqrt(nr)
			},
			restrict.loadings ={
				factors <- L%*%diag(sqr.E[1:given.d], given.d)/sqrt(nc)
				loadings  <- R*sqrt(nc)
			}
		)

	list(factors= factors, loadings= loadings, fitted.values = fitted.values
	, cov.matrix = cov.matrix, eigen.values= Eval ,Sd2 = Sd2, given.d= given.d
	, data.dim = c(nr, nc), dual=dual, L=L)		
	}


pca.fit <- function(dat, given.d=NULL 
		, restrict.mode= c("restrict.factors","restrict.loadings")
		, allow.dual = TRUE, neglect.neg.ev = TRUE){
	is.regular.panel(dat, stopper = TRUE)
	svd.pca.obj  <- svd.pca(dat, given.d = given.d, allow.dual= allow.dual
							, neglect.neg.ev = neglect.neg.ev)
	result <- restrict.pca(svd.pca.obj)
	structure(result, class = "pca.fit")
	}



