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
KSS(formula,
    consult.dim.crit = FALSE,
    additive.effects = c("none","individual","time", "twoways"),
    two.step         = FALSE,
    level            = 0.01,
    factor.dim       = NULL,
    d.max            = NULL,
    sig2.hat         = NULL,
    restrict.mode    = c("restrict.factors","restrict.loadings"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{An object of class 'formula'.}
  \item{consult.dim.crit}{logical.
    \itemize{
      \item{If \code{consult.dim.crit} is \code{FALSE} (default) and
	\code{factor.dim} is \code{NULL}:
	Only the dimensionality criterion of Kneip, Sickles & Song 2011
	is used.}
      \item{If \code{consult.dim.crit} is \code{TRUE} and
	\code{factor.dim} is \code{NULL}:
	All implemented dimensionality criteria are computed and the user
	has to select a dimension via a GUI. See also the details.}
    }
  }
  \item{additive.effects}{
    Type of Data Transformations:
    \itemize{
      \item \code{"none"}: for no transformation
      \item \code{"individual"}: for within transformation
      \item \code{"time"}: for between transformation
      \item \code{"twoways"}: for twoways transformation
  }}
  \item{two.step}{logical. If \code{two.step} is \code{FALSE} (default)
    then the original estimation procedure proposed in Kneip,
    Sickles & Song 2011 is used. If \code{two.step} is \code{TRUE} then
    the slope-parameters are re-estimated given the estimated factor
    structure. See also the details.}
  \item{level}{alpha-level for Dimensionality-Criterion of Kneip,
    Sickles & Song 2011}
  \item{factor.dim}{Dimension of Factor-Structure, pre-specified
    by the user. The default (\code{factor.dim=NULL}) yields to an internal
    estimation.}
  \item{d.max}{Maximal Dimension for Dimensionality-Estimation of Bai
    (2009). The default (\code{d.max=NULL}) yields to an internal selection of
    d.max.}
  \item{sig2.hat}{Standard deviation of the error-term, chosen by
    the user. The default (\code{sig2.hat=NULL}) yields to an internal
    estimation.}
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
      functions for a continuous argument t, as well as stochastic
      processes for discrete argument t.
      \itemize{
	\item{\code{formula}: }{Usual 'formula'-object. If you wish to
	  estimate a model without an intercept use '-1' in the
	  formula-specification. Each Variable has to be given as a
	  TxN-matrix.}
	\item{consult.dim.crit}{
	  \itemize{
	    \item{If \code{consult.dim.crit} is \code{FALSE} (default)
	      only the dimensionality criterion
	      of Kneip, Sickles & Song 2011 is used -- except for the
	      case if the user pre-specifies a dimension via the
	      argument \code{factor.dim}. If the latter is the case,
	      there is no estimation of the dimensionality of the factor
	      structure.}
	    \item{If \code{consult.dim.crit} is \code{TRUE} the user is
	      asked to select a dimension on basis of several
	      dimensionality criteria. The GUI prints the results of the
	      dimensionality criteria of Bai (???), Kneip, Sickles &
	      Song (2010), Onatski (???), and (???).}
	  }
	}
	\item{\code{additive.effects}}{
	  \itemize{
	    \item \code{"none"}: The data is not transformed, except for an
	    eventually subtraction of the overall mean, if the model is
	    estimated with an intercept. The assumed model can be written as
	    \eqn{$Y_{it}=\alpha+\sum_{j=1}^p\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}\quad i=1,\dots,n, t=1,\dots,T.$}{
              Y_{it}=\alpha+\sum_{j=1}^p\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}, i=1,...,n, t=1,...,T.}
	    The parameter 'alpha' is set to zero if '-1' is used in \code{formula}.
	    \item \code{"individual"}: This is the "within"-model
	    assuming that there are time-constant individual effects,
	    tau_i, besides the individual time trends v_i(t). The assumed model can be written as
	    \eqn{$Y_{it}=\alpha+\sum_{j=1}^p\beta_{j}
	      X_{itj}+v_i(t)+\tau_{i}+\epsilon_{it}\quad i=1,\dots,n, t=1,\dots,T.$}{
              Y_{it}=\beta_0(t)+\sum_{j=1}^p\beta_{j}
	      X_{itj}+v_i(t)+\tau_{i}+\epsilon_{it}, i=1,...,n, t=1,...,T.}
	    The parameter 'alpha' is set to zero if '-1' is used in \code{formula}.
	    \item \code{"time"}: This is the "between"-model assuming
	    that there is a common (for all individuals) time trend, beta_0(t). The
	    assumed model can be written as
	    \eqn{$Y_{it}=\alpha+\beta_{0}(t)+\sum_{j=1}^p\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}\quad i=1,\dots,n, t=1,\dots,T.$}{
              Y_{it}=\beta_0(t)+\sum_{j=1}^p\beta_{j}
	      X_{itj}+v_i(t)+\epsilon_{it}, i=1,...,n, t=1,...,T.}
	    The parameter 'alpha' is set to zero if '-1' is used in \code{formula}.
	    \item \code{"twoways"}: This is the usual twoways-model (within &
	    between) assuming that there are time-constant individual effects, tau_i, and a common (for
	    all individuals) time trend, beta_0(t). The assumed model can be written as
	    \eqn{$Y_{it}=\alpha+\beta_{0}(t)+\sum_{j=1}^p\beta_{j}
	      X_{itj}+\tau_i+v_i(t)+\epsilon_{it}\quad i=1,\dots,n, t=1,\dots,T.$}{
              Y_{it}=\beta_0(t)+\sum_{j=1}^p\beta_{j}
	      X_{itj}+\tau_i+v_i(t)+\epsilon_{it}, i=1,...,n, t=1,...,T.}
	    The parameter 'alpha' is set to zero if '-1' is used in \code{formula}.
	  }
	  \item{\code{two.step}: }{
	    \itemize{
	      \item{If \code{two.step} is \code{FALSE} (default)
		then the original estimation procedure proposed in Kneip,
		Sickles & Song 2011 is used. This procedure is allways
		consistent and shall be preferred unless the user knows
		the true dimension of the factor structure. In the
		latter case the user shall choose \code{two.step=TRUE}
		and set \code{factor.dim} equal to the known factor
		dimension.}
	      \item{If \code{two.step} is \code{TRUE} then the
		slope-parameters are re-estimated given the estimated factor
		structure. This procedure has the advantage that, if the
		selected dimension of the factor structure is true, the
		estimates of the slope-parameters are more
		efficient than those proposed in Kneip, Sickles & Song
		(2011). But, if the selected dimension of the factor
		structure is smaller than the true dimension and if the
		factor structure correlates with at least one of the
		regressors then the estimates of the slope-parameters
		suffer from an unobserved variable bias.}
	    }
	  }
	}
      }
    }
    \value{
      'KSS' returns an object of 'class' '"KSS"'.

      An object of class '"lm"' is a list containing at least the
      following components:
      \itemize{
	\item{dat.matrix: }{Whole data set stored within a
	  (N*T)x(p+1)-Matrix, where P is the number of independent
	  variables without the intercept.}
	\item{dat.dim: }{Vector of length 3: c(T,N,p)}
	\item{slope.para: }{Beta-parameters}
	\item{beta.V: }{Covariance matrix of the beta-parameters.}
	\item{names: }{Names of the dependent and independent variables.}
	\item{is.intercept: }{Used an intercept in the formula?: TRUE or FALSE}
	\item{additive.effects: }{Additive effect type. One of: "none","individual","time", "twoways".}
	\item{Intercept: }{Intercept-parameter}
	\item{Add.Ind.Eff: }{Estimated values of additive individual effects.}
	\item{Add.Tim.Eff: }{Estimated values of additive time effects.}
	\item{unob.factors: }{TxL-matrix of estimated unobserved common
	  factors, where 'L' is the number of used factors.} 
	\item{ind.loadings: }{NxL-matrix of loadings parameters.}
	\item{unob.fact.stru: }{TxN-matrix of the estimated factor
      structure. Each column represents an estimated individual unobserved time trend.}
	\item{used.dim: }{Used dimensionality of the factor structure.}
	\item{optimal.dim: }{List of proposed dimensionalities.}
	\item{fitted.values: }{Fitted values}
	\item{orig.Y: }{Original values of the dependent variable.}       
	\item{residuals: }{Residuals}    
	\item{sig2.hat: }{Estimated variance of the error term.}     
	\item{degrees.of.freedom: }{Degrees of freedom}
	\item{call}
      }
}
\references{
  \itemize{
    \item Kneip, A., Sickles, R. C., Song, W., 2011 \dQuote{A New Panel Data
      Treatment for Heterogeneity in Time Trends}, \emph{Forthcoming in
      Econometric Theory}
    \item Bai, J., 2009 \dQuote{Panel data models with interactive fixed
      effects}, \emph{Econometrica}
}}
\author{
  Dominik Liebl
}
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
KSScf, eup
}
\examples{
## See the example in 'help(Cigar)' in order to take a look at the
## data set Cigar

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
  cpi      <- matrix(Cigar$cpi, T,N)

  ## Real Disposable Income per Capita  
  l.r.ndi    <- log(matrix(Cigar$ndi,   T,N)/cpi)

  ## Real Price per Pack of Cigarettes 
  l.r.price  <- log(matrix(Cigar$price, T,N)/cpi)

  ## Real Min-Price per Pack of Cigarettes in Neighboring States
  l.r.pimin  <- log(matrix(Cigar$pimin, T,N)/cpi)

################
## ESTIMATION ##
################

## With additive time effects and 3-dimensional factor structure:
  Cigar.time <- KSS(formula          = l.sales ~ l.r.ndi + l.r.price + l.r.pimin,
                    consult.dim.crit = FALSE,
                    additive.effects = "time",
                    factor.dim       = 3)
  ## (Try also to set 'consult.dim.crit = TRUE' and 'factor.dim = NULL')

  ## Summary-Table
  (Cigar.time.summary <- summary(Cigar.time))

  ## Plot estimated common factors and individual time-trends
  plot(Cigar.time.summary)
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
