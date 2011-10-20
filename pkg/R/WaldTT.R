#Eup.object <- Eup(Y ~ X1 + X2 + X3, additive.effects = "individual")

WaldTT <- function(Eup.object, level = 0.05)
	{
	if(Eup.object$additive.effects!= "individual") stop("Additive effects should be with specified with 'individual' ")
	nr = Eup.object$dat.dim[1]
	nc = Eup.object$dat.dim[2]
	P = Eup.object$dat.dim[3]

	sig2.hat <- Eup.object$sig2.hat

	alpha <- Eup.object$Add.Ind.Eff
	alpha.bar <- mean(alpha)
	w <- sqrt(nc/2)*((nr/nc)*sum((alpha - alpha.bar)^2/sig2.hat) - 1 )
	qnor <- qnorm(1- level)
	p.value <- 1 - pnorm(w)

  if(p.value == 0) star = "***"
  else{
	if(p.value <= 0.001) star = "**"
	else {
		if(p.value <= 0.01) star = "*"
		else {
			if(p.value <= 0.05) star = "."	
			else star = " "
			}
	      }
	}
  Result <- matrix(c(round(w, 3), round(qnor, 3), (format.pval(p.value)), star), 1, 4)
  rownames(Result) <- "Existence of Individual Effects:"
  colnames(Result) <- c("w-Value", "norm-Q", "Pr( > w)", "")
  Result
	}

WaldTT(Eup(Y ~ X1 + X2 + X3, additive.effects = "individual"), level = 0.05)