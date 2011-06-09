\name{KSS}
\alias{KSS}
\alias{KSS.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
KSS-Routine
}
\description{
Estimation of Panel Data Models with Heterogeneous Time Trends
}
\usage{
KSS(formula,\dots)
\method{KSS}{default}(formula, effect = c("none", "individual", "time", "twoways"),
                      dim.crit = c("KSS.C1", "KSS.C2"), alpha = 0.01, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{An object of class "formula".}
  \item{effect}{
    Type of Data Transformations:
    \itemize{
      \item "\code{none}": for no transformation
      \item "\code{individual}": for within transformation
      \item "\code{time}": for between transformation
      \item "\code{twoways}": for twoways transformation
  }}
  \item{dim.crit}{
    Selection of the threshold for Dimensionality-Criterion:
    \itemize{
      \item "KSS.C1": (1-alpha)-quantile of the stand. normal
      distribution. This is the default.  
      \item "KSS.C2": (1-dyn.alpha)-quantile of the stand. normal
      distribution, where \eqn{$\textrm{dyn.alpha}=1-\texttt{pnorm(}\sqrt{2\,\log(\min(N,T))}\texttt{)}$}{dyn.alpha=1-pnorm(sqrt(2*log(min(N,T))))}.
    }
  }
  \item{alpha}{alpha for Dimensionality-Criterion "KSS.C1"}
  \item{\dots}{Additional arguments to be passed to low level functions.}

}
\details{
  'KSS' is a function to estimate panel data models with unobserved
  heterogeneous time trends v_i(t). The considered model in Kneip,
  Sickles & Song (2011) is given by
  \eqn{$Y_{it}=\beta_{0}(t)+\sum_{j=1}^p\beta_{j}
    X_{itj}+v_i(t)+\epsilon_{it}\quad i=1,\dots,n, t=1,\dots,T.$}{
        Y_{it}=\beta_{0}(t)+\sum_{j=1}^p\beta_{j}
	X_{itj}+v_i(t)+\epsilon_{it}, i=1,...,n, t=1,...,T.}
      Where the individual time trends, v_i(t), are assumed to come from
      a finite dimensional factor model
      \eqn{$v_i(t)=\sum_{r=1}^L\theta_{ir}g_r(t),\quad\theta_{ir}\in
	R,\quad g_r\in
	L^2[0,T].$}{v_i(t)=\sum_{r=1}^L\theta_{ir}g_r(t),
	 \theta_{ir}\in R,  g_r\in
	L^2[0,T].}
      The unobserved functions v_i(t) can be interpreted as smooth
      functions for a continuous argument t, as well as stochastic processes for discrete argument t.
      \itemize{
	\item{\code{formula}: }{Usual 'formula'-object. If you wish to
	  estimate a model without an intercept use '-1' in the
	  formula-specification. Each Variable has to be given as a TxN-matrix.}
	\item{\code{effect}}{
	  \itemize{
	    \item \code{none}: The data is not transformed, except for an
	    eventually subtraction of the overall mean, if the model is
	    estimated with an intercept. The assumed model can be written as
	    \eqn{$Y_{it}=\alpha+\sum_{j=1}^p\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}\quad i=1,\dots,n, t=1,\dots,T.$}{
              Y_{it}=\alpha+\sum_{j=1}^p\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}, i=1,...,n, t=1,...,T.}
	    The parameter 'alpha' is set to zero if '-1' is used in \code{formula}.
	    \item \code{individual}: This is the "within"-model
	    assuming that there are time-constant individual effects,
	    tau_i, besides the individual time trends v_i(t). The assumed model can be written as
	    \eqn{$Y_{it}=\alpha+\sum_{j=1}^p\beta_{j}
	      X_{itj}+v_i(t)+\tau_{i}+\epsilon_{it}\quad i=1,\dots,n, t=1,\dots,T.$}{
              Y_{it}=\beta_0(t)+\sum_{j=1}^p\beta_{j}
	      X_{itj}+v_i(t)+\tau_{i}+\epsilon_{it}, i=1,...,n, t=1,...,T.}
	    The parameter 'alpha' is set to zero if '-1' is used in \code{formula}.
	    \item \code{time}: This is the "between"-model assuming
	    that there is a common (for all individuals) time trend, beta_0(t). The
	    assumed model can be written as
	    \eqn{$Y_{it}=\alpha+\beta_{0}(t)+\sum_{j=1}^p\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}\quad i=1,\dots,n, t=1,\dots,T.$}{
              Y_{it}=\beta_0(t)+\sum_{j=1}^p\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}, i=1,...,n, t=1,...,T.}
	    The parameter 'alpha' is set to zero if '-1' is used in
	    \code{formula}.
	    \item \code{twoways}: This is the usual twoways-model (within &
	    between) assuming that there are time-constant individual effects, tau_i, and a common (for
	    all individuals) time trend, beta_0(t). The assumed model can be written as
	    \eqn{$Y_{it}=\alpha+\beta_{0}(t)+\sum_{j=1}^p\beta_{j}
	      X_{itj}+\tau_i+v_i(t)+\epsilon_{it}\quad i=1,\dots,n, t=1,\dots,T.$}{
              Y_{it}=\beta_0(t)+\sum_{j=1}^p\beta_{j}
	      X_{itj}+\tau_i+v_i(t)+\epsilon_{it}, i=1,...,n, t=1,...,T.}
	    The parameter 'alpha' is set to zero if '-1' is used in \code{formula}.
	  }
	}
      }
    }
    \value{
 \item{mu}{Intercept-Parameter as T-vector, where each element of the
   T-vector is equal to the scalar 'alpha'.}
 \item{tau}{N-vector that contains the time-constant individual effects,
   tau_i.}
 \item{beta.0}{T-vector that contains the individual-common time
   effects, beta_0(t).}
 \item{sig2.hat}{Estimated variance of the error term epsilon_{it}}
 \item{beta.V}{Var- covariance matrix of the beta-parameters.}
 \item{fitted.values}{The fitted values as TxN-matrix}
 \item{Orig.Y}{Original dependent observations, y_{ti}, as TxN-matrix}
 \item{beta}{Estimated beta-parameters}
 \item{effect}{The chosen effect-type.}
 \item{fAFactMod}{Object of the class 'fAFactMod'.}
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
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.
% }
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
