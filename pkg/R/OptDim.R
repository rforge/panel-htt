##rm(list=ls())
##################################### Bai ##################################
bai.dim.opt <- function(Obj, d.max = NULL, sig2.hat = NULL)
	{
	nr    <- Obj$nr
	nc    <- Obj$nc
	w 	<- Obj$V.d/(nr*nc)
	d.seq <- Obj$d.seq
	C     <- min(nr, nc)
	max.rk <- length(w)
	d.max <- ifelse(is.null(d.max),min(trunc(0.25*C), max.rk ), d.max)
	if(d.max > max.rk)
		{
		warning(expression("The given maximun dimension 'd.max' is larger than the postive eingen values of the corresponding covariance matrix. We use maximum rank as 'd.max' "))
		d.max <- max.rk
		}

	ref   <- ifelse(is.null(sig2.hat),w[(d.max+1)], sig2.hat)
	rmw   <- w - ref
	dmax  <- length(rmw[rmw>= 0]) - 1

	PC1 = w + ref * d.seq * ((nr+nc)/(nr*nc)) * log((nr*nc)/(nr + nc))
	d.opt.PC1 <- order(PC1)[1]-1	
	PC2 = w + ref *d.seq * ((nr+nc)/(nr*nc)) * log(C)
	d.opt.PC2 <- order(PC2)[1]-1	
	PC3 = w + ref * d.seq * (log(C)/C)
	d.opt.PC3 <- order(PC3)[1]-1	

	IC1 = log(w[1:(d.max+1)]) + d.seq[1:(d.max+1)] * ((nr+nc)/(nr*nc)) * log((nr*nc)/(nr + nc))
	d.opt.IC1 <- order(IC1)[1]-1	
	IC2 = log(w[1:(d.max+1)]) + d.seq[1:(d.max+1)] * ((nr+nc)/(nr*nc)) * log(C)
	d.opt.IC2 <- order(IC2)[1]-1	
	IC3 = log(w[1:(d.max+1)]) + d.seq[1:(d.max+1)] * (log(C)/C)
	d.opt.IC3 <- order(IC3)[1]-1	

	IPC1 = w + ref * d.seq * (nr/(4*log(log(nr)))) *((nr+nc)/(nr*nc)) * log((nr*nc)/(nr + nc))
	d.opt.IPC1 <- order(IPC1)[1]-1	
	IPC2 = w + ref *d.seq * (nr/(4*log(log(nr)))) *((nr+nc)/(nr*nc)) * log(C)
	d.opt.IPC2 <- order(IPC2)[1]-1
	IPC3 = w + ref * d.seq * (nr/(4*log(log(nr)))) *((nr+nc - ref)/(nr*nc)) *(log(nr*nc))
	d.opt.IPC3 <- order(IPC3)[1]-1

	result <- matrix(c(d.opt.PC1, d.opt.PC2, d.opt.PC3, d.opt.IC1
				, d.opt.IC2, d.opt.IC3, d.opt.IPC1, d.opt.IPC2
				, d.opt.IPC3, w[d.opt.PC1+1], w[d.opt.PC2+1]
				, w[d.opt.PC3+1], w[d.opt.IC1 + 1], w[d.opt.IC2+1]
				, w[d.opt.IC3+1], w[d.opt.IPC1+1], w[d.opt.IPC2+1]
				, w[d.opt.IPC3+1]), 9,2)
	all.crit <- c("PC1", "PC2",  "PC3", "IC1"
				, "IC2", "IC3", "IPC1", "IPC2", "IPC3")
	Result <- data.frame(I(all.crit), result
				, rep.int(ref, 9), rep(min(dmax, d.max), 9))
	colnames(Result) <- c("Criterion ", "Optimal Dimension", "sd2"
				, "sd2.hat.ref", "d.ref.max")
	return(Result)
	}


B.OptDim <- function(Obj, criteria = c("PC1","PC2","PC3","IC1","IC2","IC3"
			, "IPC1","IPC2","IPC3") , d.max = NULL, sig2.hat = NULL){
	# what is Obj?
	if(class(Obj)=="svd.pca"|class(Obj)=="fsvd.pca") obj <- Obj
	else{
		if(class(Obj)=="pca.fit"|class(Obj)=="fpca.fit"){
			nr  <- Obj$data.dim[1]
			nc  <- Obj$data.dim[2]
			V.d <- Obj$Sd2*(nr*nc)
			d.seq = seq.int(0, (length(V.d)-1))
			obj <- list(V.d = V.d, nr = nr, nc = nc, d.seq = d.seq)
			}	
		else{
			if(is.regular.panel(Obj)) obj <- svd.pca(Obj)
			else{
				if(!is.vector(Obj[[1]])|!is.numeric(Obj[[1]])
				  |!is.numeric(Obj[[2]])|length(Obj[[2]])!=2)
				  stop(c("'Obj' does not have the correct form."))
				else{# the function can deal with a list containing a vector of RSS for each d in the first listcomponent and the dimension as 2 dimensional vector in the second component 
					nr <- Obj[[2]][1]
					nc <- Obj[[2]][2]
					V.d <- Obj[[1]]
					d.seq = seq.int(0, (length(V.d)-1))
					obj <- list(V.d =V.d, nr = nr, nc = nc
						, d.seq = d.seq)
					}
				}
			}
		}
	result <- bai.dim.opt(obj, d.max = d.max, sig2.hat = sig2.hat)
	criteria <- match.arg(criteria, several.ok = TRUE)
	return(result[result[,1] %in% criteria, ])
	}
## #### Test
## Obj <- dat
## B.OptDim(pca.fit(Obj), d.max =3)

## Obj <- svd.pca(dat)
## B.OptDim(Obj, criteria = c("IPC3"), d.max = 3)

## pcaObj <- svd.pca(dat)
## Obj <- list(pcaObj$V.d, c(pcaObj$nr, pcaObj$nc))
## B.OptDim(Obj)


#####################################################################################################################

onatski.dim.opt <- function(svd.pca.obj, d.max = NULL)
	{
	nr      <- svd.pca.obj$nr
	nc      <- svd.pca.obj$nc
	w 	  <- svd.pca.obj$V.d/(nr*nc)
	ev      <- svd.pca.obj$E/nr
	max.rk  <- length(ev)

	exa.ev  <- c(sum(ev), ev)
	d.max   <- ifelse(is.null(d.max), max.rk-5, min(d.max, (max.rk-5)))	
	j       <- max(1, (d.max - 1))

	c.reg   <- as.matrix(seq((j - 1), (j + 3))^{(2/3)}, 4, 1)
	delta   <- 2* abs(coef(lm.fit(c.reg, ev[j:(j+4)])))
	dist.ev <- exa.ev[1:(max.rk-1)] -  exa.ev[2:max.rk] - delta

	FUN.o.dim <- function(o.dist.ev, d)
		{
		if(o.dist.ev[d] < 0) (d - 2)
		else {
			if (d == (max.rk-1)){
				(d - 2)
				warning(expression("ED faild to find an appropriate dimension")) 
				}
			else FUN.o.dim(o.dist.ev, (d+1))
			}
		}
	
	d.opt.i <- FUN.o.dim(dist.ev, 1)

	if(d.opt.i==d.max) {
			result <- data.frame(I("ED")
					, matrix(c(d.opt.i, w[d.opt.i+1]), 1, 2))
			colnames(result) <- c("Criterion", "Optimal Dimension", "sd2")
			return(result)
			}
	else {
		if (d.opt.i == (max.rk-1)) {
			result <- data.frame(I("ED")
					, matrix(c(d.opt.i, w[d.opt.i+1]), 1, 2))
			colnames(result) <- c("Criterion", "Optimal Dimension", "sd2")
			return(result)
			warning(expression("ED faild to find an appropriate dimension")) 
			}
		else onatski.dim.opt(svd.pca.obj, d.opt.i)
		}
	}

O.OptDim <- function(Obj, d.max = NULL){
	# what is Obj?
	if(class(Obj)=="svd.pca"|class(Obj)=="fsvd.pca") obj <- Obj
	else{
		if(class(Obj)=="pca.fit"|class(Obj)=="fpca.fit"){
			nr  <- Obj$data.dim[1]
			nc  <- Obj$data.dim[2]
			V.d <- Obj$Sd2*(nr*nc)
			E   <- Obj$eigen.values*(nr*nc)
			obj <- list(V.d = V.d, nr = nr, nc = nc, E = E)
			}	
		else{
			if(is.regular.panel(Obj)) obj <- svd.pca(Obj)
			else{
				if(!is.vector(Obj[[1]])|!is.numeric(Obj[[1]])
				  |!is.numeric(Obj[[2]])|length(Obj[[2]])!=2)
				  stop(c("'Obj' does not have the correct form."))
				else{# the function can deal with a list containing a vector of RSS for each d in the first listcomponent and the dimension as 2 d-vector 'c(nr, nc)' in the second component 
					nr  <- Obj[[2]][1]
					nc  <- Obj[[2]][2]
					V.d <- Obj[[1]][-length(Obj[[1]])]
					E   <- -diff(Obj[[1]]-Obj[[1]][1])
					obj <- list(V.d = V.d, nr = nr, nc=nc, E = E)
					}
				}
			}
		}
	result <- onatski.dim.opt(obj, d.max = d.max)
	return(result)
	}

## #### Test
## Obj <- dat
## O.OptDim(Obj)
## O.OptDim(svd.pca(Obj))
## O.OptDim(pca.fit(Obj))

## pcaObj <- svd.pca(dat)
## Obj2 <- list(pcaObj$V.d, c(pcaObj$nr, pcaObj$nc))
## O.OptDim(Obj2)


#####################################################################################################################



RH.dim.opt <- function(svd.pca.obj, d.max = NULL)
	{
	nr      <- svd.pca.obj$nr
	nc      <- svd.pca.obj$nc
	w 	  <- svd.pca.obj$V.d/(nr*nc)
	ev      <- svd.pca.obj$E/nr
	d.seq   <- svd.pca.obj$d.seq
	max.rk  <- length(w)

	exa.ev  <- c(sum(ev), ev)	
	C       <- min(nr, nc)
	d.max   <- ifelse(is.null(d.max),min(trunc(0.5*C), max.rk), d.max)
	if(d.max > max.rk)
		{
		warning(expression("The given maximun dimension 'd.max' is larger than the postive eingen values of the corresponding covariance matrix. We use maximum rank as 'd.max' "))
		d.max <- max.rk
		}

	ER = exa.ev[1:(d.max+1)]/exa.ev[2:(d.max+2)]
	d.opt.ER <- order(ER, na.last = FALSE)[d.max]
	GR = (log(w[1:(d.max+1)]) - log(w[2:(d.max+2)]))/(log(w[2:(d.max+2)]) - log(w[3:(d.max+3)]))
	d.opt.GR <- order(GR, na.last = FALSE)[d.max]

	result <- matrix(c(d.opt.ER, d.opt.GR, w[d.opt.ER+1], w[d.opt.GR+1]), 2, 2)
	Result <- data.frame(I(c("ER", "GR")),  result, rep(d.max, 2))
	colnames(Result) <- c("Criterion", "Optimal Dimension", "sd2", "d.ref.max")
	return(Result)
	}


RH.OptDim <- function(Obj, criteria = c("ER", "GR"), d.max = NULL){
	# what is Obj?
	if(class(Obj)=="svd.pca"|class(Obj)=="fsvd.pca") obj <- Obj
	else{
		if(class(Obj)=="pca.fit"|class(Obj)=="fpca.fit"){
			nr    <- Obj$data.dim[1]
			nc    <- Obj$data.dim[2]
			V.d   <- Obj$Sd2*(nr*nc)
			E     <- Obj$eigen.values*(nr*nc)
			d.seq <- seq.int(0, (length(V.d)-1))
			obj <- list(V.d = V.d, nr = nr, nc = nc, E = E)
			}	
		else{
			if(is.matrix(Obj)) obj <- svd.pca(Obj)
			else{
				if(!is.vector(Obj[[1]])|!is.numeric(Obj[[1]])
				  |!is.numeric(Obj[[2]])|length(Obj[[2]])!=2)
				  stop(c("'Obj' does not have the correct form."))
				else{# the function can deal with a list containing a vector of RSS for each d in the first listcomponent and the dimension as 2 d-vector 'c(nr, nc)' in the second component 
					nr  <- Obj[[2]][1]
					nc  <- Obj[[2]][2]
					V.d <- Obj[[1]][-length(Obj[[1]])]
					E   <- -diff(Obj[[1]]-Obj[[1]][1])
					d.seq = seq.int(0, (length(Obj[[1]])-1))
					obj <- list(V.d = V.d, nr = nr, nc = nc
						, E = E, d.seq = d.seq)
					}
				}
			}
		}

	result <- RH.dim.opt(obj, d.max = d.max)
	criteria <- match.arg(criteria, several.ok = TRUE)
	return(result[result[,1] %in% criteria, ])
	}

## #### Test
## Obj <- dat
## RH.OptDim(Obj, criteria = c("GR"))
## RH.OptDim(svd.pca(Obj))
## Obj <- pca.fit(dat)
## RH.OptDim(Obj, d.max = 10)

## pcaObj <- svd.pca(dat)
## Obj <- list(pcaObj$V.d, c(pcaObj$nr, pcaObj$nc))
## RH.OptDim(Obj)


#####################################################################################################################
KSS.dim.opt <- function(obj, sig2.hat = NULL, alpha=0.01, d.max = NULL){
  # kann direkt mit fsvd.pca-objekten arbeiten
  nr       <- obj$nr
  nc       <- obj$nc
  spar.low <- obj$spar.low
  dat      <- obj$Q.orig
  dat.smth <- obj$Q.orig.smth
  w        <- obj$V.d/(nr*nc)  
  evec     <- obj$L
  Eval     <- obj$V.d
  max.rk   <- length(Eval)
  
### calculate traces
  
  I.smth1 <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = diag(rep(1,nr)), spar = spar.low, method = 1)$ysmth
  I.smth2 <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = I.smth1,         spar = spar.low, method = 1)$ysmth
  tr.dim.zero     <- sum(diag(I.smth2))
  tr.dim.zero.sqr <- sum(diag(I.smth2)^2)

  P             <- diag(1, nr) - tcrossprod(evec)
  pca.fit.p.obj <- eigen(P)
  W             <- pca.fit.p.obj[[2]]
  ## is left out, since P.E[1:T]==rep(1,T)
  P.E           <- pca.fit.p.obj[[1]]

  W.smth  <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = W,               spar = spar.low, method = 1)$ysmth
  I.smth1 <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = diag(rep(1,nr)), spar = spar.low, method = 1)$ysmth
  I.smth2 <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = I.smth1,         spar = spar.low, method = 1)$ysmth
  tr.dim.zero     <- sum(diag(I.smth2))
  tr.dim.zero.sqr <- sum(diag(I.smth2)^2)

  diag.Wsmt <- diag(crossprod(W.smth))

  tr1 <- c(tr.dim.zero,     (sum(diag.Wsmt)   - cumsum(diag.Wsmt)))
  tr2 <- c(tr.dim.zero.sqr, (sum(diag.Wsmt^2) - cumsum(diag.Wsmt^2)))


### determine / calculate sig2.hat

  if(is.null(sig2.hat)){
	# estimation of sig2.hat: Classical, if one wants to use the KSS-Criterion for non-smoothed dat
	if(is.null(dat.smth)| !is.null(d.max)){
  		if(is.null(d.max)) d.max <- round(sqrt(min(nr, nc)))
		sig2.hat <- w[d.max+1]*(nc*nr)/(nr*nc - (nr + nc)*d.max - 1)
	}
	# estimation of sig2.hat: Variance-Estimator see Section 3.4 (KSS) 
	else{
		tr		<- (nr + sum(diag(I.smth2)) - 2*sum(diag(I.smth1)))
		sig2.hat	<- sum((dat-dat.smth)^2)/((nc-1)*tr)
                
#	rice np.variance estimator
#		sig2.hat	<- sum(diff(dat)^2)/(2*(nr - 1)*nc)	
	}

  }

## calculate the criteria of all dimensions:
  delta       <- (Eval - (nc-1) * sig2.hat * tr1[1:max.rk])/(sig2.hat * sqrt(2*nc*tr2[1:max.rk]))
  thres1      <- qnorm(1-alpha)
  thres2      <- sqrt(2*log(min(nr, nc)))# default alpha = NULL / falls alpha != 0 dann, werden beide beide berechnet
  level2      <- 1 - pnorm(thres2)
  crit1       <- delta - thres1
  crit2       <- delta - thres2
  d.opt.KSS1  <- length(crit1[crit1 > 0])# minus 1, weil start bei dim = 0
                                         # plus  1, weil nur die dim, die das crit nicht erfüllen.
  d.opt.KSS2  <- length(crit2[crit2 > 0])
  result1     <- c(d.opt.KSS1, w[d.opt.KSS1+1], sig2.hat, alpha )
  result2     <- c(d.opt.KSS2, w[d.opt.KSS2+1], sig2.hat, level2)# level2: p.value of thres2 (thres2: alternativ crit.value)

  result      <- rbind(result1, result2)
  Result      <- vector("list", 2)
  Result[[1]] <- data.frame(I(c("KSS.C1", "KSS.C2")), result)
  colnames(Result[[1]]) <- c("Criterion", "Optimal Dimension", "sd2.rest", "sd2.hat", "level")
  rownames(Result[[1]]) <- c("KSS.1", "KSS.2")
##   Result[[2]] <- c(round(delta[d.opt.KSS1],2), round(1-pnorm(delta[d.opt.KSS1]), 2), round(thres1, 2), round(alpha, 2))
##   names(Result[[2]]) <- c("Test.Stat", "p.value", "crit.value", "sig.level")
  Result[[2]] <- list(Test.Stat  = round(delta[d.opt.KSS1],2), p.value = round(1-pnorm(delta[d.opt.KSS1]), 2),
                      crit.value = round(thres1, 2), sig.level = round(alpha, 2))
return(Result)
}

## KSS.OptDim() =====================================================================================
## Called by: fAFactMod()
## Calls    : KSS.dim.opt()
##==========================

KSS.OptDim <- function(Obj,
                       criteria    = c("KSS.C1", "KSS.C2"),
                       sig2.hat    = NULL,
                       alpha       = 0.01, 
                       spar.low    = NULL,
                       d.max       = NULL){
  ## what is Obj?
  if(class(Obj)=="svd.pca"|class(Obj)=="fsvd.pca"){
    if(class(Obj)=="fsvd.pca") obj <- Obj
    else{
      ## Liste um spar.low und Q.orig.smth erweitern:
      nr          <- Obj$nr
      nc          <- Obj$nc
      spar.low    <- 0          
      Q.orig      <- Obj$Q.orig
      Q.orig.smth <- NULL
      L           <- Obj$L
      V.d         <- Obj$V.d  
      obj <- list(nr = nr, nc = nc, spar.low = spar.low, Q.orig = Q.orig,
                  Q.orig.smth = Q.orig.smth, L = L, V.d = V.d)
    }    
  }
  if(class(Obj)=="pca.fit"|class(Obj)=="fpca.fit"){
    if(class(Obj)=="fpca.fit"){
      ## umbenennungen zu (f)svd.pca-Elementen
      nr          <- Obj$data.dim[1]
      nc          <- Obj$data.dim[2]
      spar.low    <- Obj$spar.low
      Q.orig      <- Obj$orig.values
      Q.orig.smth <- Obj$orig.values.smth
      L           <- Obj$L
      V.d         <- Obj$Sd2*(nr*nc)      
      obj <- list(nr = nr, nc = nc, spar.low = spar.low, Q.orig = Q.orig,
                  Q.orig.smth = Q.orig.smth, L = L, V.d = V.d)
    }else{
      nr          <- Obj$data.dim[1]
      nc          <- Obj$data.dim[2]
      spar.low    <- 0
      Q.orig      <- Obj$orig.values
      Q.orig.smth <- NULL
      L           <- Obj$L
      V.d         <- Obj$Sd2*(nr*nc)     
      obj <- list(nr = nr, nc = nc, spar.low = spar.low, Q.orig = Q.orig,
                  Q.orig.smth = Q.orig.smth, L = L, V.d = V.d)
    } 
  }else{
    if(is.matrix(Obj)){      
      ## fPCA
      fsvd.obj  <- fsvd.pca.ramsay(Q = Obj)
    }
  }
  result        <- KSS.dim.opt(fsvd.obj, sig2.hat = sig2.hat, alpha = alpha, d.max = d.max)
  criteria      <- match.arg(criteria, several.ok = TRUE)
  Result        <- vector("list", 2)
  Result[[1]]   <- result[[1]][result[[1]][,1] %in% criteria, ]
  Result[[2]]   <- result[[2]]
  return(Result)
}
##===========================================================================================================



############## compare  ######################################################################################

OptDim.default <- function(Obj, criteria.of = c("Bai", "KSS", "Onatski", "RH")
				, d.max = NULL, sig2.hat=NULL, level= 0.05
				, spar = NULL){
	criteria.of <- match.arg(criteria.of, several.ok = TRUE)
	FUN.crit <- function(criteria.of, d.max, sig2.hat, level, spar) {
		switch(criteria.of,
		Bai = { B.OptDim(Obj = Obj, d.max = d.max, sig2.hat=NULL)
			},
		KSS = { KSS.OptDim(Obj = Obj, sig2.hat = sig2.hat, alpha = level, spar = spar)[[1]]
			},
		Onatski = {O.OptDim(Obj = Obj, d.max  = d.max)
			},
		RH = { RH.OptDim(Obj = Obj, d.max = d.max)
			})
		}

	structure(sapply(criteria.of, FUN.crit, d.max = d.max, sig2.hat=NULL
		, level= level, spar = spar, simplify = FALSE)
		, class = "OptDim")

	}

OptDim <- function(x,...) UseMethod("OptDim")

summary.OptDim <- function(x, ...)
	{
	CritNames <- c(x$Bai[,1], x$KSS[, 1], x$Onatski[,1], x$RH[, 1])
	Critdim   <- c(x$Bai[,2], x$KSS[, 2], x$Onatski[,2], x$RH[, 2])
	CritDim   <- matrix(Critdim, 1, length(Critdim))
	colnames(CritDim) <- CritNames
	rownames(CritDim) <- " "
	return(CritDim)
	}


EstDim <- function(	Obj, 
				dim.criterion = c("PC1", "PC2", "PC3",
							"IC1", "IC2", "IC3",
						      "IPC1","IPC2", "IPC3",
							"KSS.C1", "KSS.C2",
                                          "ED",  "ER",  "GR"),
				d.max,
				sig2.hat,
				level = 0.01,
				spar
			)
	{
  ## missing parameters

  	 if(missing(d.max))      d.max       <- NULL
 	 if(missing(sig2.hat))   sig2.hat    <- NULL
  	 if(missing(spar))   	 spar    <- NULL

  ## estimation 
	 dim.criterion <- match.arg(dim.criterion)
 	 est.dim       <- switch(dim.criterion,
                          PC1  = B.OptDim(Obj, criteria = c("PC1")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          PC2  = B.OptDim(Obj, criteria = c("PC2")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          PC3  = B.OptDim(Obj, criteria = c("PC3")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          IC1  = B.OptDim(Obj, criteria = c("IC1")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          IC2  = B.OptDim(Obj, criteria = c("IC2")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          IC3  = B.OptDim(Obj, criteria = c("IC3")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          IPC1 = B.OptDim(Obj, criteria = c("IPC1")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          IPC2 = B.OptDim(Obj, criteria = c("IPC2")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          IPC3 = B.OptDim(Obj, criteria = c("IPC3")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          ED   = O.OptDim(Obj, d.max = d.max),				
                          ER   = RH.OptDim(Obj, criteria = c("ER")
                            , d.max = d.max),
                          GR   = RH.OptDim(Obj, criteria = c("GR")
                            , d.max = d.max),
                          KSS.C1  = KSS.OptDim(Obj, criteria = c("KSS.C1")
				    , sig2.hat = sig2.hat, alpha=level, spar.low= spar)[[1]],
                          KSS.C2  = KSS.OptDim(Obj, criteria = c("KSS.C2")
				    , sig2.hat = sig2.hat, alpha=level, spar.low= spar)[[1]]
                          )
	est.dim
	}


## #### Test
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/Package_Version_31_3_2010/Generate_FPCAData.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel-htt/pkg/R/fpca.fit.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel-htt/pkg/R/pca.fit.R")
## ## create data for FPCA
## library(pspline)
## rm(dat)
## FS.dat      <- sim.3dim.fpca.equi(T = 300, N = 50, dim=4, sig.error = 0.07*(1/N^{0.25}), class = "matrix")

## OptDim(FS.dat, criteria.of = "KSS")

## Obj <- svd.pca(FS.dat)
## OptDim(Obj, d.max = 10)

## pcaObj <- svd.pca(dat)
## Obj <- list(pcaObj$V.d, c(pcaObj$nr, pcaObj$nc))
## OptDim(dat)
## #################################################
