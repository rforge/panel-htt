OptDim.default <- function(Obj, 
                   criteria = c("PC1", "PC2", "PC3",
                     "IC1", "IC2", "IC3",
                     "IPC1","IPC2", "IPC3",
                     "KSS.C",
                     "ED",  "ER",  "GR"),
			standardize = FALSE,
			d.max,
			sig2.hat,
			level = 0.01
                   )
  {
  ## if Eup or KSS classes
	if(class(Obj)=="KSS"| class(Obj)== "Eup") Obj <- Obj$unob.fact.stru + Obj$residuals
	else is.regular.panel(Obj, stopper = TRUE)

  ## missing parameters
  	 if(missing(d.max))      d.max       <- NULL
 	 if(missing(sig2.hat))   sig2.hat    <- NULL

  ## standardize
	 if(standardize) Obj <- standardize(Obj)

  ## estimation 
	 criteria <- match.arg(criteria, several.ok = TRUE)
         
 	 est.dim       <- function(criteria, d.max, sig2.hat, level){
				switch(criteria,
                          PC1  = try(B.OptDim(Obj, criteria = c("PC1")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          PC2  = try(B.OptDim(Obj, criteria = c("PC2")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          PC3  = try(B.OptDim(Obj, criteria = c("PC3")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IC1  = try(B.OptDim(Obj, criteria = c("IC1")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IC2  = try(B.OptDim(Obj, criteria = c("IC2")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IC3  = try(B.OptDim(Obj, criteria = c("IC3")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IPC1 = try(B.OptDim(Obj, criteria = c("IPC1")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IPC2 = try(B.OptDim(Obj, criteria = c("IPC2")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          IPC3 = try(B.OptDim(Obj, criteria = c("IPC3")
                            , d.max = d.max, sig2.hat = sig2.hat)),
                          ED   = try(O.OptDim(Obj, d.max = d.max)),				
                          ER   = try(RH.OptDim(Obj, criteria = c("ER")
                            , d.max = d.max)),
                          GR   = try(RH.OptDim(Obj, criteria = c("GR")
                            , d.max = d.max)),
                          KSS.C  = try(KSS.OptDim(Obj, criteria = c("KSS.C1")
				    , sig2.hat = sig2.hat, alpha=level
                             , factor.dim= NULL, d.max=d.max)[[1]])
                          )
	}

	Result <- sapply(criteria, est.dim, d.max = d.max, sig2.hat = sig2.hat, level = level, simplify = FALSE)
	summary <- sapply(1:length(criteria), function(l) Result[[l]][,1:2])
	Summary <- matrix(as.numeric(summary[2,]), ncol(summary))
	colnames(Summary) <- " "
	rownames(Summary) <- summary[1, ]
	Summary  <- t(Summary )
	Result$summary <- Summary 
	Result$criteria <- criteria
	Result$BaiNgC <- criteria%in%c("PC1", "PC2", "PC3", "IC1", "IC2", "IC3")
	Result$BaiC   <- criteria%in%c("IPC1", "IPC2", "IPC3")
	Result$KSSC    <- criteria%in%c("KSS.C")
	Result$OnatC  <- criteria%in%c("ED")
	Result$RHC    <- criteria%in%c("ER", "GR")
	Result$cl     <- match.call()
	Result$obj    <- Obj
	structure(Result, class = "OptDim")
}


# ####################### Methods ########################
	OptDim <- function(Obj, ...){ UseMethod("OptDim")}


## Print
print.OptDim <- function(x,...){
   	cat("Call: ")
	cl <- x$cl
   	print(cl)

	ll <- length(x$criteria)

	if(sum(x$KSSC) > 0){
	cat("\n---------")
	cat("\nCriterion of Kneip, Sickles, and Song (2009):\n\n")
		KSS2009 <- numeric(0)
		for(l in 1:ll) if(x$KSSC[l]) KSS2009 <- rbind(KSS2009, x[[l]][,1:2])
		dimKSS2009 <- matrix(KSS2009[, 2], 1, length(KSS2009[, 2]))
		colnames(dimKSS2009) <- "KSS.C"
		rownames(dimKSS2009) <- " "
		print(dimKSS2009)
	}
	
	if(sum(x$BaiNgC) > 0){
	cat("\n---------")
	if(sum(x$BaiNgC) > 1) cat("\nCriteria of Bai and Ng (2002):\n\n")
	else  cat("\nCriterion of Bai and Ng (2002):\n\n")
		BN2002 <- numeric(0)
		for(l in 1:ll) if(x$BaiNgC[l]) BN2002 <- rbind(BN2002, x[[l]][,1:2])
		dimBN2002 <- matrix(BN2002[, 2], 1, length(BN2002[, 2]))
		colnames(dimBN2002) <- BN2002[, 1]
		rownames(dimBN2002) <- " "
		print(dimBN2002)
	}
	if(sum(x$RHC) > 0){
	cat("\n--------")
	if(sum(x$RHC) > 1) cat("\nCriteria of Ahn and Horenstein (2008):\n\n")
	else cat("\nCriterion of Ahn and Horenstein (2008):\n\n") 
		AH2008 <- numeric(0)
		for(l in 1:ll) if(x$RHC[l]) AH2008 <- rbind(AH2008, x[[l]][,1:2])
		dimAH2008 <- matrix(AH2008[, 2], 1, length(AH2008[, 2]))
		colnames(dimAH2008) <- AH2008[, 1]
		rownames(dimAH2008) <- " "
		print(dimAH2008)
	}
	if(sum(x$BaiC) > 0){
	cat("\n---------")
	if(sum(x$BaiC) > 1) cat("\nCriteria of Bai (2004):\n\n") 
	else cat("\nCriterion of Bai (2004):\n\n")
		B2004 <- numeric(0)
		for(l in 1:ll) if(x$BaiC[l]) B2004 <- rbind(B2004, x[[l]][,1:2])
		dimB2004 <- matrix(B2004[, 2], 1, length(B2004[, 2]))
		colnames(dimB2004) <- B2004[, 1]
		rownames(dimB2004) <- " "
		print(dimB2004)
	}
	if(sum(x$OnatC) > 0){
	cat("\n---------")
	cat("\nCriterion of Onatski (2009):\n\n")
		O2009 <- numeric(0)
		for(l in 1:ll) if(x$OnatC[l]) O2009 <- rbind(O2009, x[[l]][,1:2])
		dimO2009 <- matrix(O2009[, 2], 1, length(O2009[, 2]))
		colnames(dimO2009) <- O2009[, 1]
		rownames(dimO2009) <- " "
		print(dimO2009)
	}

 }
	
## Print
summary.OptDim <- function(x,...){
   	cat("Call: ")
	cl <- x$cl
   	print(cl)

	ll <- length(x$criteria)

	if(sum(x$KSSC) > 0){
	cat("\n---------")
	cat("\nSequential testing of Kneip, Sickles, and Song (2009):\n\n")
		KSS2009 <- numeric(0)
		for(l in 1:ll) if(x$KSSC[l]) KSS2009 <- rbind(KSS2009, x[[l]])
		dimKSS2009 <- KSS2009[, c(2, 6)]
		colnames(dimKSS2009) <- c("Estimate", "Level")
		rownames(dimKSS2009) <- " KSS.C"
		print(round(dimKSS2009,3))
		cat("\n\nUsed Std Err: ")
		cat(round(KSS2009[1, 5], 3))
	}
	
	if(sum(x$BaiNgC) > 0){
	cat("\n\n---------")
	if(sum(x$BaiNgC) > 1) cat("\nCriteria of Bai and Ng (2002):\n\n")
	else  cat("\nCriterion of Bai and Ng (2002):\n\n")
		BN2002 <- numeric(0)
		for(l in 1:ll) if(x$BaiNgC[l]) BN2002 <- rbind(BN2002, x[[l]])
		dimBN2002 <- BN2002[, 2:3]
		colnames(dimBN2002) <- c("Estimate", "Std Err")
		rownames(dimBN2002) <- BN2002[, 1]
		print(round(dimBN2002, 3))
		cat("\n\nUsed d.max:")
		cat(BN2002[1, 5])
	}
	if(sum(x$RHC) > 0){
	cat("\n\n--------")
	if(sum(x$RHC) > 1) cat("\nCriteria of Ahn and Horenstein (2008):\n\n")
	else cat("\nCriterion of Ahn and Horenstein (2008):\n\n") 
		AH2008 <- numeric(0)
		for(l in 1:ll) if(x$RHC[l]) AH2008 <- rbind(AH2008, x[[l]])
		dimAH2008 <- AH2008[, 2:3]
		colnames(dimAH2008) <-  c("Estimate", "Std Err")
		rownames(dimAH2008) <- AH2008[,1]
		print(round(dimAH2008, 3))
		cat("\n\nUsed d.max:")
		cat(AH2008[1, 4])	
	}

	if(sum(x$BaiC) > 0){
	cat("\n\n---------")
	if(sum(x$BaiC) > 1) cat("\nCriteria of Bai (2004):\n\n") 
	else cat("\n\nCriterion of Bai (2004):\n\n")
		B2004 <- numeric(0)
		for(l in 1:ll) if(x$BaiC[l]) B2004 <- rbind(B2004, x[[l]])
		dimB2004 <- B2004[, 2:3]
		colnames(dimB2004) <- c("Estimate", "Std Err")
		rownames(dimB2004) <- B2004[, 1]
		print(round(dimB2004, 3))
		cat("\n\nUsed d.max:")
		cat(B2004[1, 5])
	}

	if(sum(x$OnatC) > 0){
	cat("\n\n---------")
	cat("\nCriterion of Onatski (2009):\n\n")
		O2009 <- numeric(0)
		for(l in 1:ll) if(x$OnatC[l]) O2009 <- rbind(O2009, x[[l]])
		dimO2009 <- O2009[, 2:3]
		colnames(dimO2009) <- c("Estimate", "Std Err")
		rownames(dimO2009) <- O2009[1, 1]
		print(round(dimO2009, 3))
	}


 }

plot.OptDim <- function(x, main, border, col, ...){
	Resultcrit <- x$summary	
	frq <- hist(Resultcrit, plot= FALSE, right = FALSE, breaks =seq.int(0, (max(Resultcrit)+1)))
	dims <- seq.int(0, max(Resultcrit))[frq$intensities != 0]


	dcol <- rainbow(length(dims))

	Obj <- x$obj
	nr <- nrow(Obj)
	nc <- ncol(Obj)
	dual = FALSE
	if(nr > nc) dual = TRUE
	if(dual)  Q <- t(Obj)%*%Obj
	else Q <- Obj%*%t(Obj)
	eig <- eigen(Q, only.values = TRUE)[[1]]
	d.max = min(nc, nr)
	perceig <- eig/sum(eig)
	perccum <- cumsum(perceig)
	ycords <- c((2*perceig[1]-perceig[2]), perceig[-d.max])
	yycoreds <- ycords 
	yycoreds[1] <- (ycords[1] - perceig[1])/2 + perceig[1]


	critperdim <- sapply(dims, function(x) paste(colnames(Resultcrit)[Resultcrit[1,]== x], collapse = ",    "))


	if(missing(main)) main = "Screeplot"
	if(missing(border)) border = FALSE
	if(missing(col)) col = "lightslategrey"
	xaxis <- barplot(perceig,  main = main, border = border 
		,axes = FALSE, col = col, ylim = c(0, ycords[1])
		, ylab = "Proportion of variance"
		, xlab = "Ordered eigenvalues")
	
	axis(1, xaxis[1:max(dims)], 1:max(dims), tick = FALSE)
	axis(2, perceig[1:max(dims)], paste((round(perceig[1:max(dims)], 3)*100), "%"))
	axis(4, yycoreds[(dims+1)], dims)

	for(r in 1:length(dims)){
	if(dims[r]!=0) abline(h = yycoreds[(dims[r]+1)], lty = 3, col = "gray")
	}
	text(rep(d.max, length(dims)), (yycoreds[(dims+1)]+ 0.01) , critperdim, adj = 1, cex = 1.1)
	}



