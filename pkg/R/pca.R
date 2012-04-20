
pcaFA.default <- function(dat, factor.dim = NULL,
	center = c("factors", "individuals", "twoways", "non"),
	standardize = FALSE,  
	consult.crit = c("Bai", "KSS", "Onatski", "RH"), 
	use.crit = c("PC1", "PC2", "PC3",
                     "IC1", "IC2", "IC3",
                     "IPC1","IPC2", "IPC3",
                     "KSS.C1", "KSS.C2",
                     "ED",  "ER",  "GR"),
	rotate = c("non", "rot.factors", "rot.loadings"),
	rot.method = varimax, eps = 1e-05,
	d.max = NULL, sig2.hat=NULL, level= 0.05, 
	restrict.mode = c("restrict.factors","restrict.loadings"), 
	allow.dual = TRUE, neglect.neg.ev = TRUE){

	is.regular.panel(dat, stopper = TRUE)

	nr <- nrow(dat)
	nc <- ncol(dat)

	center <- match.arg(center)#center = "non"
	center <- switch(center,  non = "non"
					, individuals = "individual"
					, factors = "time"
					, twoways = "twoways") 
	TransDat <- FUN.with.trans(dat, nc, nr, is.intercept = FALSE
					, effect = center)
	trDat <- TransDat$TDM
	fMean <- TransDat$TRm$TiVC
	IMean <- TransDat$TRm$InC

	if(standardize) trDat <- diag(diag(tcrossprod(trDat)/nr)^{-.5}
										, nr)%*%trDat
	
	svd.pca.obj  <- svd.pca(trDat, given.d = factor.dim
				, allow.dual= allow.dual
				, neglect.neg.ev = neglect.neg.ev)

	criteria.of <- match.arg(consult.crit, several.ok = TRUE)# criteria.of = c("Bai", "KSS", "Onatski", "RH")
	OptDimen <- OptDim(svd.pca.obj, criteria = criteria.of
			, d.max = d.max, sig2.hat = sig2.hat, level = level)

	
	if(is.null(factor.dim)){
		use.crit <- match.arg(use.crit)#use.crit ="PC1"
		factor.dim <- try(EstDim(svd.pca.obj, dim.criterion = use.crit, 
		, d.max = d.max, sig2.hat = sig2.hat, level = level)[[2]])
		if(!is.numeric(factor.dim)) factor.dim = 0
	}

	restrict.mode <-match.arg(restrict.mode)# restrict.mode = "restrict.loadings"
	rescal.result <- restrict.pca(svd.pca.obj, restrict.mode= restrict.mode)
	factors <- rescal.result$factors[, 0:factor.dim, drop = FALSE]
	loadings <- rescal.result$loadings[, 0:factor.dim, drop = FALSE]
	fitted.values <- factors%*%t(loadings) + fMean + 
						matrix(IMean, nr, nc, byrow= TRUE)
	exp.var <- svd.pca.obj$E[0:factor.dim]/sum(svd.pca.obj$E)

	rotFUN <- match.fun(rot.method)
	rotate <- match.arg(rotate)
	switch(rotate, 
			rot.factors = {
			rotmat   <- rotFUN(factors, eps = eps)$rotmat
			factors  <- factors%*%rotmat
			loadings <- loadings%*%rotmat
			exp.var  <-  t(rotmat)%*%diag(exp.var, factor.dim)%*%rotmat
			exp.var  <- diag(exp.var)
			Or <- order(exp.var, decreasing = TRUE)
			factors  <- factors[, Or]
			loadings <- loadings[, Or]
			exp.var <- exp.var[Or]
			},
			rot.loadings = {
			rotmat   <- rotFUN(loadings, eps = eps)$rotmat
			factors  <- factors%*%rotmat
			loadings <- loadings%*%rotmat
			exp.var  <-  t(rotmat)%*%diag(exp.var, factor.dim)%*%rotmat
			exp.var  <- diag(exp.var)
			Or <- order(exp.var, decreasing = TRUE)
			factors  <- factors[, Or]
			loadings <- loadings[, Or]
			exp.var <- exp.var[Or]
			},
			non = {
			factors  <- factors
			loadings <- loadings
			exp.var  <- exp.var
			}
	)
	if((factor.dim+1)> length(rescal.result$Sd2)) sd2.hat <- 0
	else sd2.hat <- rescal.result$Sd2[(factor.dim+1)]
	structure(list(
			OptDimen = OptDimen, 
			covmat = rescal.result$cov.matrix, 
			eigen.values = rescal.result$eigen.values,
			bestfactors = rescal.result$factors, 
			bestloadings = rescal.result$loadings, 
			factors= factors, 
			loadings = loadings,
			exp.var = exp.var, 
			fMean = fMean,
			IMean = IMean, 
			transformed.data = trDat,
			fitted.values = fitted.values, 
			input.data = dat,
			sd2.hat = sd2.hat, 
			factor.dim = factor.dim
		), class = "pcaFA")
}

pcaFA <- function(x, ...) UseMethod("pcaFA")

summary.pcaFA <- function(x, ...){
	OptDimen <- x$OptDimen
	cat("\nDimensiality Criteria:\n\n")
	print(summary(OptDimen))
	cat("\n\ Used Dimension       :", x$factor.dim)
	cat("\n\ Estimatd var of resid:", x$sd2.hat, "\n")
	cat("\nProportions of the explained variance:","\n\n")
	nam <- paste("Comp", 1:x$factor.dim, sep = "")
	exp.var <- x$exp.var
	prop <- rbind(exp.var, cumsum(exp.var))
	if(x$factor.dim == 0) cat(" ", 0, "\n\n")
	else{
	colnames(prop) <- nam
	rownames(prop) <- c("Relative  :", "Cumulative:")
	print(round(prop, 4))
	}
	}

plot.pcaFA <- function(x, ...){
	nr <-nrow(x$input.data)
	nc <- ncol(x$input.data)
	exp.var <- round(x$exp.var, 2)
	evv <- rep(0, min(nr, nc))
	ev <- x$eigen.values
	evv[0:length(ev)] <- ev 
	barplot(evv, axes = FALSE, main =  "Scree Plot"
					, names.arg = 1:min(nr,nc), border = NA)
	axis(4, evv[0:length(exp.var)], exp.var)
	axis(2, evv, round(evv, 2))
	mtext(expression(Variances), side = "2", line = "3")
	}

biplot.pcaFA <- function(x, choice = 1:2, sign = c(1,1), ...){
	nr <-nrow(x$input.data)
	nc <- ncol(x$input.data)
	rNam <- rownames(x$input.data)
	cNam <- colnames(x$input.data)
	if(is.null(rNam)) rNam <- paste("f", 1:nr, sep = "")
	if(is.null(cNam)) cNam <- paste("s", 1:nc, sep = "")
	factor.dim <- x$factor.dim 

	AxesNam <- c("Comp 1", "Comp 2")

	rotfactors <- x$factors[, choice]%*%diag(sign, 2)
	orifactors <- x$bestfactors[, choice]%*%diag(sign, 2)	
	
	rotloadings <- x$loadings[, choice]%*%diag(sign, 2)
	oriloadings <- x$bestloadings[, choice]%*%diag(sign, 2) 
	
	par(mfcol = c(2, 2))
	plot(rotloadings[,1], rotloadings[,2], main = "rotated scores", typ = "n"
		, xlab = AxesNam[choice[1]], ylab = AxesNam[choice[2]])
	abline(h = 0); abline(v= 0)
	text(rotloadings[,1], rotloadings[,2], cNam)
	plot(oriloadings[,1], oriloadings[,2], main = "unrotated scores", typ = "n"
		, xlab = AxesNam[choice[1]], ylab = AxesNam[choice[2]])
	abline(h = 0); abline(v= 0)
	text(oriloadings[,1], oriloadings[,2], cNam)
	plot(rotfactors[,1], rotfactors[,2], main = "rotated factors", typ = "n"
		, xlab = AxesNam[choice[1]], ylab = AxesNam[choice[2]])
	abline(h = 0); abline(v= 0)
	text(rotfactors[,1], rotfactors[,2], rNam) 
	plot(orifactors[,1], orifactors[,2], main = "unrotated factors", typ = "n"
		, xlab = AxesNam[choice[1]], ylab = AxesNam[choice[2]])
	abline(h = 0); abline(v= 0)
	text(orifactors[,1], orifactors[,2], rNam) 
	}
