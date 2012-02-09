\name{KSScf}
\alias{KSScf}
\alias{KSScf.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
KSS Common-Factors Specificationtest
}
\description{
Test for certain user-specified common factors
}
\usage{
KSScf(formula,
      g.function,
      additive.effects = c("none", "individual", "time", "twoways"),
      level)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{An object of class "formula".}
  \item{g.function}{User-specified common factors given as a vector,
    matrix, or list of functions.}
    \item{additive.effects}{
    Type of Data Transformations:
    \itemize{
      \item \code{"none"}: for no transformation
      \item \code{"individual"}: for within transformation
      \item \code{"time"}: for between transformation
      \item \code{"twoways"}: for twoways transformation
  }}
  \item{level}{alpha-level for Dimensionality-Criterion of Kneip,
    Sickles & Song 2011}
  }
  \details{
    'KSScf' is a function to test for certain user-specified common
    factors. The user can give his common factors to the argument
    \code{g.function}.
    \itemize{
      \item{\code{g.function}: }{Can be a vector of length T, a matrix of
	dimension TxL, or a list of length L with executable functions
	in each list-element.}
    }
}
\references{
  \itemize{
    \item Kneip, A., Sickles, R. C., Song, W., 2011 \dQuote{A New Panel Data
      Treatment for Heterogneity in Time Trends}, \emph{Forthcoming in
      Econometric Theory}
}}
\author{
  Dominik Liebl
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

## Unobserved common factors:
CF.1 <- DGP1$CF.1[1:T]
CF.2 <- DGP1$CF.2[1:T]
CF.3 <- DGP1$CF.3[1:T]

#############################################
## Test for user-specified common factors: ##
#############################################

## H0: 'The given factors are the true common factors'

  ## Cannot be rejected:
  KSScf(formula    = Y ~ X1 + X2,
        g.function = matrix(c(CF.1,CF.2,CF.3), nrow=T,ncol=3))

  ##  Can be rejected:
  KSScf(formula    = Y ~ X1 + X2,
        g.function = CF.1)
}

