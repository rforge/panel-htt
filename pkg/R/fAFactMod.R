######################################################################

fAFactMod <- function(dat, dim.criterion = c("KSS", "PC1", "PC2", "PC3",  "IC1",
                                             "IC2", "IC3", "IPC1","IPC2", "IPC3",
                                             "ED",  "ER",  "GR"),
                      factor.dim , d.max,
                      alpha=0.01,
                      sig2.hat,
                      sig2.hat.mode = c("functional", "classical"),  
                      restrict.mode = c("restrict.factors","restrict.loadings"), 
                      allow.dual = TRUE)
{
  
  nr   <- nrow(dat)
  nc   <- ncol(dat)
  
  # missing parameters

  if(missing(factor.dim)) factor.dim  <- NULL
  if(missing(d.max))      d.max       <- NULL
  if(missing(sig2.hat))   sig2.hat    <- NULL

  # smoothing the residuals (small degree of undersmoothing)

  spar.low <- smooth.Pspline(x=seq(0, 1, length.out=nr), y=dat, method = 3       )$spar * 0.8
  dat.smth <- smooth.Pspline(x=seq(0, 1, length.out=nr), y=dat, spar   = spar.low)$ysmth

  # fpca.fit with smoothed data

  fpca.fit.obj <- fpca.fit(dat           = dat.smth,
                           given.d       = factor.dim,
                           restrict.mode = restrict.mode,
                           allow.dual    = allow.dual)
        
  # estimation of sig2.hat:      

  # Variance-Estimator Section 3.4 (KSS)
  if(dim.criterion=="KSS" & sig2.hat.mode=="functional"){
    id.smth1  <- smooth.Pspline(x = seq.int(0,1, length.out= nr) , y = diag(1,nr),  spar = spar.low)$ysmth
    id.smth2  <- smooth.Pspline(x = seq.int(0,1, length.out= nr) , y = id.smth1,    spar = spar.low)$ysmth
    tr        <- (nr + sum(diag(id.smth2)) - 2*sum(diag(id.smth1)))
    sig2.hat  <- sum((dat-dat.smth)^2)/((nc-1)*tr)
  }
  
  # Variance-Estimator traditional
  if(dim.criterion=="KSS" & sig2.hat.mode=="classical"){
    # Note: no smoothing of the residuals for var-estimation
    sig2.hat <- svd.pca(dat)$V.d[1]/(nr*nc)
  }
        
  # dimension selection
  dim.criterion <- match.arg(dim.criterion)
  est.dim       <- switch(dim.criterion,
                          PC1  = B.OptDim(fpca.fit.obj, criteria = c("PC1")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          PC2  = B.OptDim(fpca.fit.obj, criteria = c("PC2")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          PC3  = B.OptDim(fpca.fit.obj, criteria = c("PC3")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          IC1  = B.OptDim(fpca.fit.obj, criteria = c("IC1")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          IC2  = B.OptDim(fpca.fit.obj, criteria = c("IC2")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          IC3  = B.OptDim(fpca.fit.obj, criteria = c("IC3")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          IPC1 = B.OptDim(fpca.fit.obj, criteria = c("IPC1")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          IPC2 = B.OptDim(fpca.fit.obj, criteria = c("IPC2")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          IPC3 = B.OptDim(fpca.fit.obj, criteria = c("IPC3")
                            , d.max = d.max, sig2.hat = sig2.hat),
                          ED   = O.OptDim(fpca.fit.obj, d.max = d.max),				
                          ER   = RH.OptDim(fpca.fit.obj, criteria = c("ER")
                            , d.max = d.max),
                          GR   = RH.OptDim(fpca.fit.obj, criteria = c("GR")
                            , d.max = d.max),
                          KSS  = KSS.OptDim(fpca.fit.obj, sig2.hat = sig2.hat,
                            alpha=alpha,  spar.low = spar.low)
                          )
  print(est.dim)
  opt.d  <- est.dim[1,2]
  used.d <- ifelse(is.null(factor.dim), opt.d, factor.dim)

  # factors and loadings parameters

  if(used.d!=0){
    factors  <- fpca.fit.obj$factors[,  1:used.d, drop= FALSE]
    loadings <- fpca.fit.obj$loadings[, 1:used.d, drop= FALSE]
    dat.fit  <- tcrossprod(factors, loadings)
    sd2      <- fpca.fit.obj$Sd2[used.d+1]
    
    result        <- list(fitted.values = dat.fit,
                          factors       = factors, 
                          loadings      = loadings,
                          sd2           = sd2,
                          given.fdim    = factor.dim,
                          optimal.fdim  = opt.d,
                          used.fdim     = used.d)
  }
  else{
    factors <- matrix(0, nr, 1)
    scores  <- matrix(0, nc, 1)
    dat.fit <- tcrossprod(factors, scores)
    sd2     <- fpca.fit.obj$Sd2[used.d+1]
    
    result       <- list(fitted.values = dat.fit,
                         factors       = factors,
                         loadings      = loadings,
                         resid.sd2     = sd2,
                         given.fdim    = factor.dim, 
                         optimal.fdim  = opt.d,
                         used.fdim     = used.d)
  }
  return(result)
}
