\name{Spatial Weights}
\docType{data}
\alias{W}
\title{Spatial Weights for Cigarette Consumtion Data}
\description{
  A matrix of dimension 46x46, holding the spatial weights for the 46
  States in the Data set of Baltagi, B.H. and  D. Levin (1992). 
}

\usage{data(W)}

\format{A data frame containing :
  \describe{
    \item{W}{spatial weights for the 46 US States of the Cigar Data}
  }
}

\references{
Baltagi, B.H. and  D. Levin (1992) \dQuote{Cigarette taxation: Raising revenues and reducing consumption}, \emph{Structural Changes and Economic Dynamics}, \bold{3}, 321--335.
}

\keyword{datasets}
\examples{
data(W)
W.mat <- matrix(W, nrow=46, ncol=46)
W.mat
}
