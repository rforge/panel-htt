\name{checkSpecif}
\alias{checkSpecif}
\alias{checkSpecification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tests to check the Model Specifications
}
\description{
Tests to check the Model Specifications
}
\usage{
checkSpecif(obj1, obj2, level = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj1}{Bla}
  \item{obj2}{Bla}
  \item{level}{Bla}
}
%\details{}
%\value{}
\references{
  \itemize{
    \item Bai, J., 2009 \dQuote{Panel data models with interactive fixed
      effects}, \emph{Econometrica}
    \item Kneip, A., Sickles, R. C., Song, W., 2012 \dQuote{A New Panel Data
      Treatment for Heterogeneity in Time Trends}, \emph{Econometric Theory}
}}
\author{
  Oualid Bada, Dominik Liebl
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
