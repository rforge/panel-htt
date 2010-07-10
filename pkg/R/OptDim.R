
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
KSS.dim.opt <- function(obj, sig2.hat = NULL, alpha=alpha){
  nr      <- obj$nr
  nc      <- obj$nc
  w       <- obj$V.d/(nr*nc)  
  evec    <- obj$Evec
  Eval    <- obj$V.d
  max.rk  <- length(Eval)

  ## estimation of sig2.hat:      
  if(is.null(sig2.hat)){
    ## Variance-Estimator see Section 3.4 (KSS) =============================================================
    id.smth1  <- smooth.Pspline(x = seq.int(0,1, length.out= nr) , y = diag(1,nr),  spar = spar.low)$ysmth
    id.smth2  <- smooth.Pspline(x = seq.int(0,1, length.out= nr) , y = id.smth1,    spar = spar.low)$ysmth
    tr        <- (nr + sum(diag(id.smth2)) - 2*sum(diag(id.smth1)))
    sig2.hat  <- sum((dat-dat.smth)^2)/((nc-1)*tr)
    ##=======================================================================================================
  }
P             <- diag(1, nr) - tcrossprod(evec)
pca.fit.p.obj <- eigen(P)
W             <- pca.fit.p.obj[[2]]
P.E           <- pca.fit.p.obj[[1]]

W.smth  <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = W,               spar = spar.low, method = 1)$ysmth
I.smth1 <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = diag(rep(1,nr)), spar = spar.low, method = 1)$ysmth
I.smth2 <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = I.smth1,         spar = spar.low, method = 1)$ysmth
tr.dim.zero     <- sum(diag(I.smth2))
tr.dim.zero.sqr <- sum(diag(I.smth2)^2)

diag.Wsmt <- diag(crossprod(W.smth))

tr1 <- c(tr.dim.zero,     (sum(diag.Wsmt)   - cumsum(diag.Wsmt)))
tr2 <- c(tr.dim.zero.sqr, (sum(diag.Wsmt^2) - cumsum(diag.Wsmt^2)))

delta     <- (Eval - (nc-1) * sig2.hat * tr1[1:max.rk])/(sig2.hat * sqrt(2*nc*tr2[1:max.rk]))
thres     <- qnorm(1-alpha)
crit      <- delta - thres
d.opt.KSS <- length(crit[crit > 0])# minus 1, weil start bei dim = 0
                                   # plus 1, weil nur die dim, die das crit nicht erfüllen.

  result <- data.frame(I("KSS"), matrix(c(d.opt.KSS, w[d.opt.KSS+1]), 1, 2))
  colnames(result) <- c("Criterion", "Optimal Dimension", "sd2")
  print(result)
  return(result)
}



KSS.OptDim <- function(Obj, sig2.hat=NULL, alpha = 0.01){
	# what is Obj?
  if(class(Obj)=="svd.pca"|class(Obj)=="fsvd.pca"){
    nr    <- Obj$nr
    nc    <- Obj$nr
    V.d   <- Obj$V.d
    E     <- Obj$E
    Evec  <- Obj$L
    
    obj <- list(V.d = V.d, nr = nr, nc = nc, E = E, Evec = Evec)
  }
  else{
    if(class(Obj)=="pca.fit"|class(Obj)=="fpca.fit"){
      nr    <- Obj$data.dim[1]
      nc    <- Obj$data.dim[2]
      V.d   <- Obj$Sd2*(nr*nc)
      E     <- Obj$eigen.values*(nr*nc)
      Evec  <- Obj$L
      d.seq <- seq.int(0, (length(V.d)-1))
      
      obj   <- list(V.d = V.d, nr = nr, nc = nc, E = E, Evec = Evec)
    }	
## 		else{
## 			if(is.matrix(Obj)) obj <- fsvd.pca(Obj)
## 			else{
## 				if(!is.vector(Obj[[1]])|!is.numeric(Obj[[1]])
## 				  |!is.numeric(Obj[[2]])|length(Obj[[2]])!=2)
## 				  stop(c("'Obj' does not have the correct form."))
## 				else{# the function can deal with a list containing a vector of RSS for each d in the first listcomponent and the dimension as 2 d-vector 'c(nr, nc)' in the second component 
## 					nr  <- Obj[[2]][1]
## 					nc  <- Obj[[2]][2]
## 					V.d <- Obj[[1]][-length(Obj[[1]])]
## 					E   <- -diff(Obj[[1]]-Obj[[1]][1])
## 					d.seq = seq.int(0, (length(Obj[[1]])-1))
## 					obj <- list(V.d = V.d, nr = nr, nc = nc
## 						, E = E, d.seq = d.seq)
## 					}
## 				}
## 			}
  }
  
  result   <- KSS.dim.opt(obj, sig2.hat, alpha, spar.low)
  return(result)
}

#############################################################################################################


############## compare

OptDim <- function(Obj, criteria.of = c("Bai", "Onatski","KSS", "RH")
				, d.max = NULL, sig2.hat=NULL, alpha= 0.05
				, spar.low = 0.005){
	criteria.of <- match.arg(criteria.of, several.ok = TRUE)
	FUN.crit <- function(criteria.of, d.max, sig2.hat, alpha, spar) {
		switch(criteria.of,
		Bai = { B.OptDim(Obj = Obj, d.max = d.max, sig2.hat=NULL)
			},
		Onatski = {O.OptDim(Obj = Obj, d.max  = d.max)
			},
		KSS = { KSS.OptDim(Obj = Obj, sig2.hat = sig2.hat, alpha = alpha, spar.low = spar.low)
			},
		RH = { RH.OptDim(Obj = Obj, d.max = d.max)
			})
		}

	structure(sapply(criteria.of, FUN.crit, d.max = d.max, sig2.hat=NULL
		, alpha= 0.05, spar = 0.005, simplify = FALSE)
		, class = "OptDim")

	}

## #### Test
## Obj <- dat
## OptDim(Obj)

## Obj <- svd.pca(dat)
## OptDim(Obj, d.max = 10)

## pcaObj <- svd.pca(dat)
## Obj <- list(pcaObj$V.d, c(pcaObj$nr, pcaObj$nc))
## OptDim(Obj)

