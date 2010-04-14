##########################
# restrict mode "standard": F'F/T = I and S'S/N = diag
# restrict mode "inverse" : S'S/N = I and F'F/T = I

restrict.pca <- function(pca.fit.obj, restric.mode= c("standard","inverse.standard"))
	{
	if(class(pca.fit.obj)!="pca.fit") stop(c("The pca.fit.obj is not a 'pca.fit' object"))
 # pca.fit object  
	nr    <- nrow(pca.fit.obj$Q.fit)
	nc    <- ncol(pca.fit.obj$Q.fit)
	U 	<- pca.fit.obj$U
	W 	<- pca.fit.obj$W
	sqr.E <- pca.fit.obj$sqr.E
	E	<- pca.fit.obj$E
	Sd2   <- pca.fit.obj$V.d/(nr*nc)

 # restrict eigen values
	Eval	<- E/(nr*nc)	

 # restric factors and scores
	re.mo <-match.arg(restric.mode)
	
	switch(re.mo,
			standard=
			{
			factors <- U*sqrt(nr)
			scores  <- W%*%diag(sqr.E)/sqrt(nr)
			},
			inverse.standard=
			{
			factors <- U%*%diag(sqr.E)/sqrt(nc)
			scores  <- W*sqrt(nc)
			}
		)

	list(factors= factors, scores= scores, Eval= Eval, Sd2 = Sd2)		
	}
