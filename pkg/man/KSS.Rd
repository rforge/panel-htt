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
KSS(x,\dots)
\method{KSS}{default}(formula, effect = c("none", "individual", "time", "twoways"),
                      dim.crit = c("KSS.C1", "KSS.C2"), alpha = 0.01, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{What the hell ? }
  \item{formula}{An object of class "formula".}
  \item{effect}{
    Type of Data Transformations:
    \itemize{
      \item "none": for no transformation
      \item "individual": for within transformation
      \item "time": for between transformation
      \item "twoways": for twoways transformation
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
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
% \value{
% %%  ~Describe the value returned
% %%  If it is a LIST, use
% %%  \item{comp1 }{Description of 'comp1'}
% %%  \item{comp2 }{Description of 'comp2'}
% %% ...
% }
\references{
  \itemize{
    \item Kneip, A., Sickles, R. C., Wonho, S., 2011 \dQuote{A New Panel Data
      Treatment for Heterogneity in Time Trends}, \emph{Forthcoming in
      Econometric Theory}
}}
% \author{
% %%  ~~who you are~~
% }
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.
% }
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
