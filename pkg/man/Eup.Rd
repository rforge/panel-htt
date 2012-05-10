\name{Eup}
\alias{Eup}
\alias{Eup.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Eup-Routine
}
\description{
Estimation of Panel Data Models with ...
}
\usage{
Eup(formula,
    additive.effects = c("none", "individual", "time", "twoways"),
    dim.criterion    = c("PC1", "PC2", "PC3", "IC1", "IC2" , "IC3",
"IPC1", "IPC2", "IPC3" , "ED"),
    d.max            = NULL,
    sig2.hat         = NULL,
    factor.dim       = NULL,
    double.iteration = TRUE,
    start.beta       = NULL,
    max.iteration    = 500,
    convergence      = 1e-6,
    restrict.mode    = c("restrict.factors", "restrict.loadings"),
    ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{An object of class 'formula'.}
  \item{additive.effects}{
    Type of Data Transformations:
    \itemize{
      \item \code{"none"}: for no transformation
      \item \code{"individual"}: for within transformation
      \item \code{"time"}: for between transformation
      \item \code{"twoways"}: for twoways transformation
  }}
  \item{dim.criterion}{Bla}
  \item{d.max}{Maximal Dimension for Dimensionality-Estimation of Bai
    (2009). The default (\code{d.max=NULL}) yields to an internal selection of
    d.max.}
  \item{sig2.hat}{Standard deviation of the error-term, chosen by
    the user. The default (\code{sig2.hat=NULL}) yields to an internal
    estimation.}
  \item{factor.dim}{Dimension of Factor-Structure, pre-specified
    by the user. The default (\code{factor.dim=NULL}) yields to an internal
    estimation.}
  \item{double.iteration}{Bla}
  \item{start.beta}{Bla}
  \item{max.iteration}{Bla}
  \item{convergence}{Bla}
  \item{restrict.mode}{
    Type of Restriction on the Factor-Structure:
    \itemize{
      \item \code{"restrict.factors"}: Factors are restricted to have an
      euclidean norm of 1.
      \item \code{"restrict.loadings"}: Factor-Loadings are restricted
      to have an euclidean norm of 1.
    }
  }
  \item{...}{Additional arguments to be passed to the low level functions.}
}
%\details{}
\value{
      'Eup' returns an object of 'class' '"Eup"' containing the
      following components:
      \itemize{
    \item{dat.matrix}
    \item{formula}
    \item{dat.dim}
    \item{slope.para}
    \item{names}
    \item{is.intercept}
    \item{additive.effects}
    \item{Intercept} 
    \item{Add.Ind.Eff} 
    \item{Add.Tim.Eff} 
    \item{unob.factors} 
    \item{ind.loadings} 
    \item{unob.fact.stru} 
    \item{used.dim}
    \item{proposed.dim}
    \item{optimal.dim}
    \item{factor.dim} 
    \item{d.max} 
    \item{dim.criterion} 
    \item{OvMeans}
    \item{ColMean} 
    \item{RowMean}
    \item{max.iteration} 
    \item{convergence}
    \item{start.beta} 
    \item{Nbr.iteration}  		
    \item{fitted.values} 
    \item{orig.Y} 
    \item{residuals} 
    \item{sig2.hat.dim}
    \item{sig2.hat} 
    \item{degrees.of.freedom}
    \item{call}
      }
}
\references{
  \itemize{
    \item Bai, J., 2009 \dQuote{Panel data models with interactive fixed
      effects}, \emph{Econometrica}
}}
\author{
  Oualid Bada
}
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
KSS
}
\examples{
## See the example in 'help(Cigar)' in order to take a look at the
## data set 'Cigar'

##########
## DATA ##
##########

data(Cigar)
N <- 46
T <- 30

## Dependent variable:

  ## Cigarette-Sales per Capita
  l.sales    <- log(matrix(Cigar$sales, T,N))

## Independent variables:

  ## Consumer Price Index
  cpi        <- matrix(Cigar$cpi, T,N)

  ## Real Disposable Income per Capita  
  l.r.ndi    <- log(matrix(Cigar$ndi,   T,N)/cpi)

  ## Real Price per Pack of Cigarettes 
  l.r.price  <- log(matrix(Cigar$price, T,N)/cpi)

  ## Real Min-Price per Pack of Cigarettes in Neighboring States
  l.r.pimin  <- log(matrix(Cigar$pimin, T,N)/cpi) 



}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
