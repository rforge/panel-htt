#Int.model <- Eup(Y ~ -1 + X1 + X2 + X3, factor.dim= 2)
#Add.model <- Eup(Y ~ -1 + X1 + X2 + X3, factor.dim= 0, additive.effects = "individual")
#summary(Int.model)

HausmanTT <- function(Int.model, Add.model, level = 0.05)
{
  inf.Int.mod <- Eup.inference(Int.model)
  beta.Int.mod <- coef(Int.model)$Slope.Coef
  inv.ZZ <- inf.Int.mod$inv.ZZ


  inf.Add.mod <- Eup.inference(Add.model)
  beta.Add.mod <- coef(Add.model)$Slope.Coef
  inv.XX <- inf.Add.mod$inv.ZZ

  Delta <- inv.ZZ - inv.XX
  sig2.hat <- inf.Int.mod$sig2.hat
  inv.Delta <- solve(Delta*sig2.hat)
  nr <- Int.model$dat.dim[1]
  nc <- Int.model$dat.dim[2]
  P  <- Int.model$dat.dim[3]
  J <- nr*nc*t(beta.Int.mod - beta.Add.mod)%*%inv.Delta%*%(beta.Int.mod - beta.Add.mod)
  qchi2 <- qchisq(1 - level, P)
  p.value <- 1 - pchisq(J, P)
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
  Result <- matrix(c(round(J, 3), round(qchi2, 3), (format.pval(p.value)), star), 1, 4)
  rownames(Result) <- "Additive vs. Interactive Effects:"
  colnames(Result) <- c("H-Value", "Chi2-Q", "Pr( > H)", "")
  Result
}

HausmanTT(Int.model, Add.model)