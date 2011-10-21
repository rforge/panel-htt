\name{Cigar}
\docType{data}
\alias{Cigar}
\title{Cigarette Consumption }
\description{
  a panel of N=46 observations each with time-dimension T=30 from 1963 to 1992

 \emph{total number of observations} :  1380

 \emph{observation} :  regional

 \emph{country} :  United States
}

\usage{data(Cigar)}

\format{A data frame containing :
  \describe{
    \item{state}{state abbreviation}
    \item{year}{the year}
    \item{price}{price per pack of cigarettes}
    \item{pop}{population}
    \item{pop16}{population above the age of 16}
    \item{cpi}{consumer price index (1983=100)}
    \item{ndi}{per capita disposable income}
    \item{sales}{cigarette sales in packs per capita}
    \item{pimin}{minimum price in adjoining states per pack of cigarettes}
  }
}

\source{
Online complements to Baltagi (2001).
 \url{http://www.wiley.com/legacy/wileychi/baltagi/}.
}

\references{
 Baltagi, Badi H. (2001) \emph{Econometric Analysis of Panel Data}, 2nd ed., John Wiley and Sons. 

 Baltagi, B.H. and  D. Levin (1992) \dQuote{Cigarette taxation: raising revenues and reducing consumption}, \emph{Structural Changes and Economic Dynamics}, \bold{3}, 321--335.

 Baltagi, B.H., J.M. Griffin and  W. Xiong (2000) \dQuote{To pool or not to pool: homogeneous versus heterogeneous estimators applied to cigarette demand}, \emph{Review of Economics and Statistics}, \bold{82}, 117--126.
}

\keyword{datasets}
\examples{
data(Cigar)
N <- 46
T <- 30

## Dependent variable:
  ## Cigarette-Sales per Capita
sales    <- matrix(Cigar$sales, T,N)     

## Independent variables:
  ## Consumer Price Index
  cpi      <- matrix(Cigar$cpi, T,N)

  ## Real Disposable Income per Capita  
  r.ndi    <- matrix(Cigar$ndi,   T,N)/cpi

  ## Real Price per Pack of Cigarettes 
  r.price  <- matrix(Cigar$price, T,N)/cpi

  ## Real Min-Price per Pack of Cigarettes in Neighboring States
  r.pimin  <- matrix(Cigar$pimin, T,N)/cpi 

########################################
## Take a look at the data set Cigar: ##
########################################
par(mfrow=c(2,2))
## Dependent variable
matplot(sales,main="Cigaret-Sales per Capita",type="l",xlab="Time",ylab="")
## Independent variables
matplot(r.ndi,main="Real Disposable Income per Capita",type="l",xlab="Time",ylab="")
matplot(r.price,main="Real Price of Cigarets per Pack",type="l",xlab="Time",ylab="")
matplot(r.pimin,main="Real Minimum-Price of Cigarets\n per Pack in Neighbouring States",type="l",xlab="Time",ylab="")
par(mfrow=c(1,1))
}
