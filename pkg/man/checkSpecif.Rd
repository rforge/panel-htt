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
  \item{obj1}{If 'obj2' is leaved unspecified, 'obj1' must be an object of class "`KSS"' or "`Eup"' where the model is estimated with a non-zero factor dimension. Otherwise, 'obj1' must be an object of class "`Eup"' where the underling model is estimated with "individual", "time", or "twoways"-effects and a factor dimension equal to zero.}
  \item{obj2}{Can be specified with an object of class "`Eup"' where the underling model is estimated without additive effects and a non-zero factor dimension.}
  \item{level}{The significance level used for the tests.}
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
KSS, Eup, OptDim
}
\examples{
## See the example in 'help(Cigar)' in order to take a look at the
## data set 'Cigar'

##########
## DATA ##
##########

data(Cigar)
## Panel-Dimensions:
N <- 46
T <- 30
## Dependent variable:
  ## Cigarette-Sales per Capita
  l.Consumption    <- log(matrix(Cigar$sales, T,N))
## Independent variables:
  ## Consumer Price Index
  cpi        <- matrix(Cigar$cpi, T,N)
  ## Real Price per Pack of Cigarettes 
  l.Price  <- log(matrix(Cigar$price, T,N)/cpi)
  ## Real Disposable Income per Capita  
  l.Income    <- log(matrix(Cigar$ndi,   T,N)/cpi)
  
## Testing the Sufficiency of Classical Additive Effects
	twoways.obj <- Eup(d.l.Consumption ~ -1 + d.l.Price + d.l.Income, factor.dim = 0, additive.effects = "twoways")
	not.twoways.obj <- Eup(d.l.Consumption ~ -1 + d.l.Price + d.l.Income, factor.dim = 7,
additive.effects = "none")
	checkSpecif(obj1 = twoways.obj, obj2 = not.twoways.obj, level = 0.01)
	
## 
	Eup.obj <- Eup(d.l.Consumption ~ -1 + d.l.Price + d.l.Income,
additive.effects = "twoways")
	checkSpecif(Eup.obj, level = 0.01)

}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
