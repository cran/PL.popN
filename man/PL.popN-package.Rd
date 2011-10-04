\name{PL.popN-package} 

\alias{PL.popN-package}

\docType{package}

\title{Closed population size estimation from capture-recapture data using a partial 
likelihood with covariates.
}

\description{\code{PL.popN} gives closed population size estimates and 
standard errors for capture-recapture partial likelihood models when capture 
propabilities depend on covariates. Models recognized by \code{PL.popN} must 
come in the form of the following logistic regression models: Generalized Linear
Models (GLM), Generalized Additive Models (GAM), Local Polynomial Models (LP), 
Generalized Linear Mixed Model (GLMM) or Simulation Extrapolation (SIMEX) models. 
The number of capture occasions (\code{tau}) and sample size (\code{D}) must 
also be specified.
}            

\details{
\tabular{ll}{
Package: \tab PL.popN\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2011-10-4\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}   

\references{Stoklosa, Hwang, Wu and Huggins. (2011). Heterogeneous 
capture-recapture models with covariates: A partial likelihood approach for 
closed populations. \emph{Biometrics}.
}

\author{Jakub Stoklosa, Wenhan Hwang and Richard Huggins.}

\keyword{package} 
