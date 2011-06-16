\name{Simulated Data}
\docType{data}
\alias{kss_twoways_with_Intercept}
\alias{kss_twoways_without_Intercept}
\alias{kss_time_with_Intercept}
\alias{kss_time_without_Intercept}
\alias{kss_individual_with_Intercept}
\alias{kss_individual_without_Intercept}
\alias{kss_none_with_Intercept}
\alias{kss_none_without_Intercept}

\title{Simulated Panel-Data Sets with Factor Structure}
\description{
  A Panel-Data Set with:
  
  \emph{time-index} : t=1,\ldots,T=100 and

  \emph{individual} : index i=1,\ldots,N=50

  \emph{total number of observations} :  500
}

\usage{
data(kss_twoways_with_Intercept)
data(kss_twoways_without_Intercept)
data(kss_time_with_Intercept)
data(kss_time_without_Intercept)
data(kss_individual_with_Intercept)
data(kss_individual_without_Intercept)
data(kss_none_with_Intercept)
data(kss_none_without_Intercept)
}

\format{A data frame containing :
  \describe{
    \item{Y}{dependent variable}
    \item{X1}{first regressor simulated as \eqn{$X_{1,it}\sim N(0,\sqrt{10})$}{X_{1,it}\sim N(0,\sqrt{10})}}
    \item{X2}{second regressor}
  }
}

\details{
  \itemize{
  \item{}
  }
}
% \source{
% }

% \references{
%  Baltagi, Badi H. (2001) \emph{Econometric Analysis of Panel Data}, 2nd ed., John Wiley and Sons. 

%  Baltagi, B.H. and  D. Levin (1992) \dQuote{Cigarette taxation: raising revenues and reducing consumption}, \emph{Structural Changes and Economic Dynamics}, \bold{3}, 321--335.

%  Baltagi, B.H., J.M. Griffin and  W. Xiong (2000) \dQuote{To pool or not to pool: homogeneous versus heterogeneous estimators applied to cigarette demand}, \emph{Review of Economics and Statistics}, \bold{82}, 117--126.
% }

\keyword{datasets}
