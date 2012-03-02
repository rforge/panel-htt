\name{ADDvsINT}
\alias{ADDvsINT}
\alias{ADDvsINT.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specificationtest
}
\description{
Test for ...
}
\usage{
ADDvsINT(obj,
         H0.ADD.effects = c("none", "individual", "time", "twoways"),
         level)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An object either of class "KSS" or of class "Eup".}
  \item{H0.ADD.effects}{
    Type of Null-Hypothesis:
    \itemize{
      \item \code{"none"}: H0 There are none additive effects
      \item \code{"individual"}: H0 There are only additive individual-effects
      \item \code{"time"}: H0 There are only additive time-effects
      \item \code{"twoways"}: H0 There are only additive twoways-effects
  }}
  \item{level}{alpha-level for Dimensionality-Criterion of Kneip,
    Sickles & Song 2011}
  }
  \details{
    'ADDvsINT' ... 
  }
\references{
  \itemize{
    \item Kneip, A., Sickles, R. C., Song, W., 2011 \dQuote{A New Panel Data
      Treatment for Heterogneity in Time Trends}, \emph{Forthcoming in
      Econometric Theory}
}}
\author{
  Oualid Bada and Dominik Liebl
}
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
fAFactMod  
}
\examples{
## See the example in 'help(DGP1)' in order to take a look at the
## simulated data set DGP1

data(DGP1)
N <- 60
T <- 30

## Observed variables:
Y  <- matrix(DGP1$Y,  nrow=T, ncol=N)
X1 <- matrix(DGP1$X1, nrow=T, ncol=N)
X2 <- matrix(DGP1$X2, nrow=T, ncol=N)

## Unobserved additive individual-effects:
add.ind      <- sample(c(1:100),N)
add.ind      <- add.ind-mean(add.ind)
add.ind      <- matrix(rep(add.ind,each=T),T,N)

## Panel-Data-Generation: 
Y            <-  add.ind + 5 * X1 - 5 * X2 + matrix(rnorm(T*N),T,N)

## KSS-Estimation:
KSS.obj <- KSS.default(formula=Y ~ X1 + X2)

#############################################
## Test                                    ##
#############################################

## Cannot be rejected
ADDvsINT(obj=KSS.obj, H0.ADD.effects="individual")
## Can be rejected
ADDvsINT(obj=KSS.obj, H0.ADD.effects="time")
}

