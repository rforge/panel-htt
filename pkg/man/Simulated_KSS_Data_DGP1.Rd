\name{Simulated Data for the KSS-Model: DGP1}
\docType{data}
\alias{DGP1}

\title{Simulated Panel-Data Set with Factor Structure}
\description{
  A Panel-Data Sets with:
  
  \emph{time-index} : t=1,\ldots,T=30

  \emph{individual-index} : index i=1,\ldots,N=100

}

\usage{
data(DGP1)
}

\format{A data frame containing :
  \describe{
    \item{Y}{dependent variable as N*T-vector}
    \item{X1}{first regressor as N*T-vector}
    \item{X2}{second regressor as N*T-vector}
    
    The time-index t is running "faster" than the individual-index i
    such that e.g. Y_{it} is ordered as:
    Y_{11},Y_{12},\ldots,Y_{1T},Y_{21},Y_{22},\ldots,Y_{NT}
  }
}

\details{
  The Panel-Data Set DPG1 is simulated according to the model:
  \eqn{$Y_{it}=\beta_{1}X_{it1}+\beta_{2}X_{it2}+v_i(t)+\epsilon_{it}\quad
    i=1,\dots,n;\quad
    t=1,\dots,T$}{Y_{it}=\beta_{1}X_{it1}+\beta_{2}X_{it2}+v_i(t)+\epsilon_{it},
    i=1,\dots,n, t=1,\dots,T}

  With time varying individual effect that follow second order polynomials:
  \eqn{$v_i(t)=\theta_{i0}+\theta_{i1}\frac{t}{T}+\theta_{i2}\left(\frac{t}{T}\right)^2$}{v_i(t)=theta_{i0}+theta_{i1}*frac{t}{T}+theta_{i2}*(frac{t}{T})^2}.

  The Regressors X_{it}=(X_{it1},X_{it2})' are simulated from a bivariate VAR model:
  \eqn{$X_{it}=R X_{i,t-1}+\eta_{it}\quad\textrm{with}\quad\eta_{it}\sim
    N(0,I_2)$}{X_{it}=R X_{i,t-1}+eta_{it} with eta_{it}~N(0,I_2)}

  In this Panel-Data Set the regressors are exogenous.
}

\references{
  \itemize{
    \item Kneip, A., Sickles, R. C., Song, W., 2011 \dQuote{A New Panel Data
      Treatment for Heterogneity in Time Trends}, \emph{Forthcoming in
      Econometric Theory}
  }
}
\author{
  Dominik Liebl
}

\keyword{datasets}
