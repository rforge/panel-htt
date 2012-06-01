\name{Eup}
\alias{Eup}
\alias{Eup.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Eup-Routine
}
\description{
Estimation of Panel Data Models with Interactive Fixed Effects. 
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
  \item{formula}{An object of class 'formula' where the arguments are matrices. The number of rows has to be equal to the temporal dimension and the number of columns has to be equal to the number of individuals. The details of model specification are given under ‘Details’.}
  \item{additive.effects}{
    Type of Data Transformations:
    \itemize{
      \item \code{"none"}: for no transformation
      \item \code{"individual"}: for within transformation
      \item \code{"time"}: for between transformation
      \item \code{"twoways"}: for twoways transformation
  }}
  \item{dim.criterion}{The dimensionality criterion to be used if \code{factor.dim} is left unspecified. The default criterion is \code{"PC1"}}.
  \item{d.max}{Maximal Dimension for Dimensionality-Estimation of Bai
    (2009). The default (\code{d.max=NULL}) yields to an internal selection of
    d.max.}
  \item{sig2.hat}{Standard deviation of the error-term required for the computation of some dimensionality criteria. The user can specify it in instead of \code{d.max}. The default (\code{sig2.hat=NULL}) yields to an internal estimation.}
  \item{factor.dim}{Dimension of Factor-Structure, pre-specified
    by the user. The default (\code{factor.dim=NULL}) yields to an internal estimation.}
  \item{double.iteration}{logical. If \code{FALSE} the update of the factor dimension \code{d} will be done simultaneously with remaining model parameters without alternating between inner and outer iteration. This may speed up computations, but the convergence is less stable than in the default setting.}
  \item{start.beta}{allows the user to give a vector of starting values for the slope parameters.}
  \item{max.iteration}{controls the maximum number of iterations. The default is '500'.}
  \item{convergence}{Convergence condition of the estimators. The default is '1e-6'.}
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
\details{
'Eup' is a function to estimate equidistant panel data models with unobserved
  multiple time varying individual effects. The considered model is given by
  \eqn{$Y_{it}= \sum_{j=1}^P\beta_{j}
    X_{itj}+ v_{it} + \epsilon_{it}\quad i=1,\dots,n; t=1,\dots,T.$}{
        Y_{it}= \sum_{j=1}^P\beta_{j}X_{itj}+v_{it}+ \epsilon_{it}, i=1,...,n; t=1,...,T.}
  Where the individual time trends, v_{it}, are assumed to come from a finite dimensional factor model: 
   \eqn{$v_{it}=\sum_{l=1}^d\lambda_{il}f_{lt},\quad \lambda_{il}, f_{lt} \in	R.$}{v_{it}=\sum_{l=1}^d\lambda_{il}f_{lt}.}
         \itemize{
	\item{\code{formula} }{Usual 'formula'-object. If you wish to
	  estimate a model without an intercept use '-1' in the
	  formula-specification. Each Variable has to be given as a
	  TxN-matrix. Missing values are not allowed.}
	\item{\code{additive.effects}}{
	  \itemize{
	    \item \code{"none"}: The data is not transformed, except for a subtraction of the overall mean, if the model is estimated with an intercept. The assumed model can be written as
	    \eqn{$Y_{it}= \mu + \sum_{j=1}^P\beta_{j}
	      X_{itj}+v_{it}+\epsilon_{it}\quad i=1,\dots,n; t=1,\dots,T.$}{
              Y_{it}=\mu+\sum_{j=1}^P\beta_{j}
	      X_{itj}+v_{it}+\epsilon_{it}, i=1,...,n; t=1,...,T.}
	    The parameter 'mu' is set to zero if '-1' is used in \code{formula}.
	    \item \code{"individual"}: This is the "within"-model, which 
	    assumes that there are time-constant individual effects,
	    alpha_i, besides the individual time trends v_{it}. The
	    model can be written as
	    \eqn{$Y_{it}= \mu +\alpha_{i} + \sum_{j=1}^P\beta_{j}
	      X_{itj}+v _{it} + \epsilon_{it} \quad i=1,\dots,n; t=1,\dots,T.$}{
              Y_{it}=\mu + \alpha_{i}+ \sum_{j=1}^P\beta_{j}
	      X_{itj} + v_{it} + \epsilon_{it}, i=1,...,n; t=1,...,T.}
	    The parameter 'mu' is set to zero if '-1' is used in \code{formula}.
	    \item \code{"time"}: This is the "between"-model, which assumes
	    that there is a common time trend (for all individuals), theta_t. The
	    model can be written as
	    \eqn{$Y_{it}=\mu+ \theta_t + \sum_{j=1}^P\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}\quad i=1,\dots,n; t=1,\dots,T.$}{
              Y_{it}=\mu+ \theta_t + \sum_{j=1}^P\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}, i=1,...,n; t=1,...,T.}
	    The parameter 'mu' is set to zero if '-1' is used in \code{formula}.
	    \item \code{"twoways"}: This is the "twoways"-model ("within" &
	    "between"), which assumes that there are time-constant
	    individual effects, alpha_i, and a common time trend,
	    theta_t. The model can be written as
	    \eqn{$Y_{it}=\mu+ \alpha_{i} + \theta_t +\sum_{j=1}^P\beta_{j}
	      X_{itj}+\tau_i+v_i(t)+\epsilon_{it}\quad i=1,\dots,n; t=1,\dots,T.$}{
              Y_{it}=\mu+ \alpha_{i} + \theta_t +\sum_{j=1}^P\beta_{j}
	      X_{itj}+\tau_i+v_i(t)+\epsilon_{it}, i=1,...,n; t=1,...,T.}
	    The parameter 'mu' is set to zero if '-1' is used in \code{formula}.
	  }}}
}

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
    \item Bada, O. and Kneip, A., 2012  \dQuote{Panel Data Models with Unobserved Multiple
Time-Varying Effects to Estimate Risk Premium of Corporate Bonds}, \emph{Working paper}, \emph{University of Bonn}
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
