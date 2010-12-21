######################################################################
rm(list=ls())
fAFactMod <- function(dat, demean   = TRUE,
                      add.effects   = c("none", "individual", "time", "twoways"),
                      dim.criterion = c("KSS.C1", "KSS.C2",
                        "PC1", "PC2", "PC3", "IC1", "IC2", "IC3", "IPC1", "IPC2", "IPC3", "ED", "ER", "GR"),
                      factor.dim, d.max, sig2.hat, spar.low,
                      alpha         = 0.01,
                      restrict.mode = c("restrict.factors","restrict.loadings"), allow.dual = TRUE)  
{
  ## checks and preparations
  is.regular.panel(dat, stopper = TRUE) 
  nr         <- nrow(dat)
  nc         <- ncol(dat)

  ## data-transformation
  with.trans <- match.arg(add.effects)
  dat.trans  <- FUN.with.trans(dat, N = nc, T = nr, is.intercept = demean, effect = with.trans)
  dat        <- dat.trans$TDM

  ## missing parameters
  if(missing(factor.dim)) factor.dim  <- NULL
  if(missing(d.max))      d.max       <- NULL
  if(missing(sig2.hat))   sig2.hat    <- NULL
  if(missing(spar.low))   spar.low    <- NULL
  
  ## fpca.fit (within fpca.fit(): "dat" is (under)smoothed)

  fpca.fit.obj <- fpca.fit(dat           = dat,
                           spar.low      = spar.low,    # if NULL: 0.8*spar.opt (GCV)
                           given.d       = factor.dim,  # if NULL: max.rk <- ifelse(neglect.neg.ev,nbr.pos.ev, length(Eval))
                           restrict.mode = restrict.mode,
                           allow.dual    = allow.dual) 

  
  ## dimension selection
  dim.criterion <- match.arg(dim.criterion)
  est.dim       <- EstDim(Obj           = dat,
                          dim.criterion = dim.criterion,
                          d.max         = d.max,
                          sig2.hat      = sig2.hat,
                          level         = alpha, spar= 3)
  
  opt.d         <- est.dim[1,2]
  used.d        <- ifelse(is.null(factor.dim), opt.d, factor.dim)

  ## factors and loadings parameters
  if(used.d!=0){
    factors  <- fpca.fit.obj$factors[,  1:used.d, drop= FALSE]
    loadings <- fpca.fit.obj$loadings[, 1:used.d, drop= FALSE]
    dat.fit  <- tcrossprod(factors, loadings)
    sd2      <- fpca.fit.obj$Sd2[used.d+1]
    
    result   <- list(fitted.values = dat.fit,
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
    
    result  <- list(fitted.values = dat.fit,
                    factors       = factors,
                    loadings      = loadings,
                    resid.sd2     = sd2,
                    given.fdim    = factor.dim, 
                    optimal.fdim  = opt.d,
                    used.fdim     = used.d)
  }
  return(result)
}

## ## TEST: =========================================================================================================

## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel-htt/pkg/R/OptDim.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel-htt/pkg/R/pca.fit.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel-htt/pkg/R/fpca.fit.R")
## source("/home/dom/Dokumente/Uni/Promotion/Panel_HTT/our_package/panel-htt/pkg/R/FUN.with.trans.R")
## source("/home/dom/Dokumente/Uni/Promotion/myRoutines/Generate_FS.R")

## ## create data for FPCA
## library(pspline)
## T   = 100
## N   = 50
## dim = 4

## ## FS-Structure
## FS.obj   <- sim.FS(T = T, N = N, dim=dim, Factors= "sin", AR =c(0,0), ar.sd = 0.25)
## FS.obs   <- FS.obj[[1]]

## # create data for FPCA
## library(pspline)
## dat        <- FS.obs

## ## OptDim.obj <- OptDim(dat, criteria.of="KSS")
## ## OptDim.obj

## fAF.obj    <- fAFactMod(dat, dim.criterion="KSS.C1")
## fAF.obj$used.fdim
