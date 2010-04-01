select.nf <-
function(obj, criteria = c("Bai", "Onatski","KSS", "RH")
    , d.max = NULL, sig2.hat, alpha= 0.05, spar = 0.005){
if(class(obj)!="pca.fit") obj <- pca.fit(obj)
criteria <- match.arg(criteria, several.ok = TRUE)
FUN.crit <- function(criteria, d.max) {
 switch(criteria,
 Bai = { bai.dim.opt(pca.fit.obj = obj, d.max = d.max)
  },
 Onatski = {onatski.dim.opt(pca.fit.obj = obj, d.max  = d.max)
  },
 KSS = { kss.dim.opt(pca.fit.obj = obj, sig2.hat = sig2.hat, alpha = alpha, spar = spar)
  },
 RH = { RH.dim.opt(pca.fit.obj = obj, d.max = d.max)
  })
 }

structure(sapply(criteria, FUN.crit, d.max = d.max, simplify = FALSE), class = "select.nf")
}

