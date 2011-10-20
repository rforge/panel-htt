Eup.object <- Eup(Y ~ X1 + X2 + X3, factor.dim = 1, additive.effects = "twoways")

TestIE <- function(Eup.object, level = 0.05)
	{
	if(Eup.object$additive.effects!= "twoways") stop("Additive effects should be with specified with 'twoways' ")
	nr = Eup.object$dat.dim[1]
	nc = Eup.object$dat.dim[2]
	P = Eup.object$dat.dim[3]
	sig2.hat <- Eup.object$sig2.hat

	dat.matrix <- Eup.object$dat.matrix
	y 	<- dat.matrix[, 1, drop = FALSE]
	x 	<- dat.matrix[,-1, drop = FALSE]

	df <- (nr-1)*(nc -1) - P
	resiHO2 <- sum((lm.fit(x, y)$residuals)^2)
	sig2h0 <- resiHO2/(df)
	Std.Err <- sig2.hat*sqrt(2/df)

	D <- (resiHO2 - df*sig2.hat)/(sqrt(2*sig2.hat^2*df))
	p.wert <- 1 - pnorm(D)
	inf.result <- cbind(sig2h0, Std.Err, D, p.wert)
	colnames(inf.result) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
	rownames(inf.result) <- ""
	printCoefmat(inf.result)
	}

TestIE(Eup(Y ~ X1 + X2 + X3, factor.dim = 2, additive.effects = "twoways"))