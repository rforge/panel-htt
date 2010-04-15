gafsa <-
  function(dat, x=NULL, given.d, dim.crit = c("PC1"
                                   , "PC2"  , "PC3" , "IC1"
                                   , "IC2"  , "IC3" , "IPC1"
                                   , "IPC2" , "IPC3", "ED"
                                   , "KSS.C", "ER"  , "GR")                  
   # arg for FUN.dim
           , d.max
           , alpha = 0.05
   # args for smooth.Pspline
           , smooth.dat = FALSE
           , spar = 0
           , method=3, ... )
  {
  # dat dimension

    nr <- nrow(dat)
    nc <- ncol(dat)

  # x axis 

    if(is.null(x)) x <- seq.int(1, nr)

  # smooth or not 

    d.crit <- match.arg(dim.crit)

    if(smooth.dat | d.crit=="KSS.C"){
      SMth <-smooth.Pspline(x=seq.int(0,1, length.out = nr), y=dat, method = method, spar = spar)
      smth.dat <- SMth$ysmth
      df <- SMth$df
      spar <- SMth$spar}
 if(d.crit=="KSS.C"){
   res.smth <- dat - smth.dat
   sig2.hat <- sum(diag(crossprod(res.smth)))/(nr*nc)
   #sig2.hat <- (1/((nc-1)*(nr-df)^2))*sum((dat-smth.dat)^2)
 }
    
    
  # pca fitting 
    if(smooth.dat| d.crit=="KSS.C") Q <- smth.dat else Q <- dat
 
    obj.pca.fit <- pca.fit(Q)

  # dimension test
    if(missing(d.max)) d.max <- NULL

    d.opt <- switch(d.crit,
                    PC1 = select.nf(obj.pca.fit, criteria = c("Bai")
                      , d.max = d.max)$Bai[1,1],
                    PC2 = select.nf(obj.pca.fit, criteria = c("Bai")
                      , d.max = d.max)$Bai[2,1],
                    PC3 = select.nf(obj.pca.fit, criteria = c("Bai")
                      , d.max = d.max)$Bai[3,1],
                    IC1 = select.nf(obj.pca.fit, criteria = c("Bai")
                      , d.max = d.max)$Bai[4,1],
                    IC2 = select.nf(obj.pca.fit, criteria = c("Bai")
                      , d.max = d.max)$Bai[5,1],
                    IC3 = select.nf(obj.pca.fit, criteria = c("Bai")
                      , d.max = d.max)$Bai[6,1],
                    IPC1 = select.nf(obj.pca.fit, criteria = c("Bai")
                      , d.max = d.max)$Bai[7,1],
                    IPC2 = select.nf(obj.pca.fit, criteria = c("Bai")
                      , d.max = d.max)$Bai[8,1],
                    IPC3 = select.nf(obj.pca.fit, criteria = c("Bai")
                      , d.max = d.max)$Bai[9,1],
                    ED   = select.nf(obj.pca.fit, criteria = c("Onatski")
                      , d.max = d.max)$Onatski[1,1],
                    KSS.C = select.nf(obj.pca.fit, criteria = c("KSS")
                      , sig2.hat = sig2.hat, alpha = alpha, spar = spar)$KSS[1,1],
                    ER   = select.nf(obj.pca.fit, criteria = c("RH")
                      , d.max = d.max)$RH[1,1],
                    GR   = select.nf(obj.pca.fit, criteria = c("RH")
                      , d.max = d.max)$RH[2,1]
                    )
    
    if(missing(given.d)) used.d <- d.opt else used.d <-  given.d
    
  # scores and factors according to the restriction and used dimension

    obj.rst.pca  <- restrict.pca(obj.pca.fit)
    all.factors  <- obj.rst.pca$factors
    all.scores   <- obj.rst.pca$scores
    all.sd2.resi <- obj.rst.pca$Sd2
    
    if(used.d!=0)
      {        
        factors <- all.factors[, 1:used.d, drop= FALSE]
        scores  <- all.scores[, 1:used.d, drop= FALSE]
        dat.fit <- tcrossprod(factors, scores)
        R <- list(dat.fit, x, factors, scores, all.sd2.resi[used.d+1], used.d, d.opt)
      }
    else
      {
        factors <- matrix(0, nr, 1)
        scores  <- matrix(0, nc, 1)
        dat.fit <- tcrossprod(factors, scores)
        R <- list(dat.fit, x, factors, scores, all.sd2.resi[used.d+1], used.d, d.opt)
      }
    R
  }

