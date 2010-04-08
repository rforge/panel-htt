#####################################################################################################################
# dim.opt-Functions:
# bai.dim.opt, onatski.dim.opt, AH.dim.opt and kss.dim.opt
#####################################################################################################################

bai.dim.opt <- function(pca.fit.obj, d.max = NULL)
	{
	nr    <- pca.fit.obj$nr
	nc    <- pca.fit.obj$nc
	w 	<- pca.fit.obj$V.d/(nr*nc)
	d.seq <- pca.fit.obj$d.seq
	C     <- min(nr, nc)
	max.rk <- length(w)
	d.max <- ifelse(is.null(d.max),min(trunc(0.25*C), max.rk ), d.max)
	if(d.max > max.rk)
		{
		warning(c("The given maximun dimension is 'd.max' is larger than the postive eingen values of the corresponding covariance matrix. We use maximum rank as 'd.max' "))
		d.max <- max.rk
		}

	ref   <- w[(d.max+1)] # ifelse(is.null(sig2.ref),w[(d.max+1)], sig2.ref )


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
	rownames(result) <- c("PC1", "PC2",  "PC3", "IC1", "IC2", "IC3", "IPC1", "IPC2", "IPC3")
	colnames(result) <- c("Optimal Dimension", "sd2")
	return(result)
	}

#####################################################################################################################
onatski.dim.opt <- function(pca.fit.obj, d.max = NULL)
	{
	nr      <- pca.fit.obj$nr
	nc      <- pca.fit.obj$nc
	w 	  <- pca.fit.obj$V.d/(nr*nc)
	ev      <- pca.fit.obj$E/nr
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
				warning(c("ED faild to find an appropriate dimension")) 
				}
			else FUN.o.dim(o.dist.ev, (d+1))
			}
		}
	
	d.opt.i <- FUN.o.dim(dist.ev, 1)

	if(d.opt.i==d.max) {
			result <- matrix(c(d.opt.i, w[d.opt.i+1]), 1, 2)
			rownames(result) <- c("ED")
			colnames(result) <- c("Optimal Dimension", "sd2")
			return(result)
			}
	else {
		if (d.opt.i == (max.rk-1)) {
			result <- matrix(c(d.opt.i, w[d.opt.i+1]), 1, 2)
			rownames(result) <- "ED"
			colnames(result) <- c("Optimal Dimension", "sd2")
			return(result)
			warning(c("ED faild to find an appropriate dimension")) 
			}
		else onatski.dim.opt(pca.fit.obj, d.opt.i)
		}
	}

#####################################################################################################################
RH.dim.opt <- function(pca.fit.obj, d.max = NULL)
	{
	nr      <- pca.fit.obj$nr
	nc      <- pca.fit.obj$nc
	w 	<- pca.fit.obj$V.d/(nr*nc)
	ev      <- pca.fit.obj$E/nr
	d.seq   <- pca.fit.obj$d.seq
	max.rk  <- length(w)

	exa.ev  <- c(sum(ev), ev)	
	C       <- min(nr, nc)
	d.max   <- ifelse(is.null(d.max),min(trunc(0.5*C), max.rk), d.max)
	if(d.max > max.rk)
		{
		warning(c("The given maximun dimension is 'd.max' is larger than the postive eingen values of the corresponding covariance matrix. We use maximum rank as 'd.max' "))
		d.max <- max.rk
		}

	ER = exa.ev[1:(d.max+1)]/exa.ev[2:(d.max+2)]
	d.opt.ER <- order(ER, na.last = FALSE)[d.max]
	GR = (log(w[1:(d.max+1)]) - log(w[2:(d.max+2)]))/(log(w[2:(d.max+2)]) - log(w[3:(d.max+3)]))
	d.opt.GR <- order(GR, na.last = FALSE)[d.max]

	result <- matrix(c(d.opt.ER, d.opt.GR, w[d.opt.ER+1], w[d.opt.GR+1]), 2, 2)
	rownames(result) <- c("ER", "GR") 
	colnames(result) <- c("Optimal Dimension", "sd2")
	return(result)
	}

#####################################################################################################################
kss.dim.opt <- function(pca.fit.obj, sig2.hat, alpha, spar = spar){
nr      <- pca.fit.obj$nr
nc      <- pca.fit.obj$nc
w 	<- pca.fit.obj$V.d/(nr*nc)  
evec    <- pca.fit.obj$U
Eval    <- pca.fit.obj$V.d#[-1]
max.rk  <- length(Eval)

P             <- diag(1, nr) - tcrossprod(evec)
pca.fit.p.obj <- eigen(P)
W             <- pca.fit.p.obj[[2]]
P.E           <- pca.fit.p.obj[[1]]

W.smth  <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = W,               spar = spar, method = 1)$ysmth
I.smth1 <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = diag(rep(1,nr)), spar = spar, method = 1)$ysmth
I.smth2 <- smooth.Pspline(x=seq.int(0,1, length.out = nr), y = I.smth1,         spar = spar, method = 1)$ysmth
tr.dim.zero     <- sum(diag(I.smth2))
tr.dim.zero.sqr <- sum(diag(I.smth2)^2)

diag.Wsmt <- diag(crossprod(W.smth))

tr1 <- c(tr.dim.zero,     (sum(diag.Wsmt)   - cumsum(diag.Wsmt)))
tr2 <- c(tr.dim.zero.sqr, (sum(diag.Wsmt^2) - cumsum(diag.Wsmt^2)))

delta <- (Eval - (nc-1)*sig2.hat*tr1[1:max.rk])/(sig2.hat*sqrt(2*nc*tr2[1:max.rk]))
thres <- qnorm(0.99999)
crit <- delta - thres
d.opt.KSS <- length(crit[crit > 0])# minus 1, weil start bei dim = 0
                                   # plus 1, weil nur die dim, die das crit nicht erfüllen.

result <- matrix(c(d.opt.KSS, w[d.opt.KSS+1]), 1, 2)
	  rownames(result) <- c("KSS") 
	  colnames(result) <- c("Optimal Dimension", "sd2")
	  return(result)
}

#####################################################################################################################
#kss.dim.opt1 < - function(pca.fit.obj, alpha=0.00001, sig2.hat=NULL)
#                        
#{
#	nr     <- pca.fit.obj$nr     # T
#	nc     <- pca.fit.obj$nc     # N
#	w      <- pca.fit.obj$V.d    # upper sums of eigenvalues for each l: sum(lambda_l+1,...,lambda_d.max). Ohne Umskalierung 
#	d.seq  <- pca.fit.obj$d.seq
#	max.rk <- length(d.seq)
# 	
#	sig2.hat <- ifelse(is.null(sig2.hat), (w[5]/(nc*nr)), sig2.hat) 
#	sig2.hat <- w[5]/(nc*nr)
#	tr.1     <- nr - d.seq
#	tr.2     <- tr.1
#	delta    <- (w - (nc - 1) * tr.1* sig2.hat )/(sig2.hat * sqrt(2 * nc * tr.2))
#	thres    <- log(min(nr, nc))#qnorm((1-alpha))
#	crit     <- delta - thres
#	d.opt	   <- length(crit[crit>0])
#d.opt
#}



