\name{OptDim}
\alias{OptDim}
\alias{OptDim.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimation of the Factor Dimension
}
\description{
  Functions for the Estimation of the Factor Dimension
}
\usage{
OptDim(Obj, 
       criteria = c("PC1", "PC2", "PC3",
                     "IC1", "IC2", "IC3",
                     "IPC1","IPC2", "IPC3",
                     "KSS.C",
                     "ED",  "ER",  "GR"),
       standardize = FALSE,
       d.max,
       sig2.hat,
       spar,
       level = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Obj}{The function requires either a Txn matrix or an object with class "`Eup"' or "`KSS"'.}
  \item{criteria}{A character vector that contains the desired criteria to be used. If it is left unspecified, the function returns the result of all 13 criteria. }
  \item{standardize}{logical. If \code{TRUE} the input variable will be standardized. Default is \code{FALSE}.}
  \item{d.max}{Maximal dimension used in the dimensionality-criteria of Bai
    (2009). The default (\code{d.max=NULL}) yields to an internal selection of
    d.max.}
  \item{sig2.hat}{The squared standard deviation of the error-term required for the computation of some dimensionality criteria. The user can specify it in instead of \code{d.max}. The default (\code{sig2.hat=NULL}) yields to an internal estimation.}
  \item{spar}{Smoothing parameter used to calculate the criterion of Kneip, Sickles, and Song (2012). The default is \code{NULL}, which leads to internal computation.}
  \item{level}{The significance level used for the criterion of Kneip, Sickles, and Song (2012). The default is 0.01.}
}
\details{The function 'OptDim' allows for a comparison of the optimal factor dimensions obtained from different panel criteria (in total 13) adjusted for the case where N and T are }
\value{
'OptDim' returns an object of 'class' '"OptDim"' containing a list with the
      following components:
\item{summary:} 
\item{criteria:}
\item{BaiNgC}
\item{BaiC}
\item{KSSC}
\item{OnatC}
\item{RHC}
\item{cl}
\item{obj}

}
\references{
  \itemize{
      \item  Ahn, S. C., Horenstein, A. R. 2009 \qQuote{Eigenvalue ratio test for the number of factors}, \emph{Working Paper}
    \item Bai, J., 2009 \dQuote{Panel data models with interactive fixed
      effects}, \emph{Econometrica}
    \item Bai, J. 2004 \qQuote{Estimating cross-section common stochastic trends in nonstationary data},  \emph{Journal of Econometrics}
    \item Bai, J., Ng, S. 2009 \qQuote{Determining the number of factors in approximated factor models}, \emph{Econometrica}
    \item Kneip, A., Sickles, R. C., Song, W., 2012 \dQuote{A New Panel Data
      Treatment for Heterogeneity in Time Trends}, \emph{Forthcoming in
      Econometric Theory}
    \item Onatski, A. 2010 \qQuote{Determining the number of factors from empirical distribution of eigenvalues}, \emph{The Review of Economics and Statistics}
}}
\author{
  Oualid Bada
}
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{plot.OptDim, 
% }
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
