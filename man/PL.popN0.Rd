\name{PL.popN0}

\alias{PL.popN0}

\title{Closed population size estimation from capture-recapture data using a 
partial likelihood with covariates.}

\description{\code{PL.popN} gives closed population size estimates and 
standard errors for capture-recapture partial likelihood models when capture 
propabilities depend on covariates. Models recognized by \code{PL.popN} must 
come in the form of the following logistic regression models: Generalized Linear
Models (GLM), Generalized Additive Models (GAM), Local Polynomial Models (LP), 
Generalized Linear Mixed Model (GLMM) or Simulation Extrapolation (SIMEX) models. 
The number of capture occasions (\code{tau}) and sample size (\code{D}) must 
also be specified.
}

\usage{PL.popN0(model, tau, D)}

\arguments{
\item{model}{ an object of model type, i.e. model objects from either 
\code{glm, gam, sm.binomial, glmmPQL, simex}. }
\item{tau}{ the number of capture occasions in the experiment. }
\item{D}{ the number of unique individuals captured in the experiment. }
}

\details{The structure of the partial likelihood allows the user to construct 
logistic regression models and analyse the data when considering: parametric 
models when using \code{glm()}, non-parametric (semi-parametric) smoothing when 
using  \code{gam()} (from the \code{mgcv} package) or \code{sm.binomial()}, 
mixed models by using \code{glmmPQL} and measurement error when usng 
\code{simex}. Population size estimates are obtained by extracting fitted 
capture probabilities from the models, which are then employed in the 
Horwitz-Thompson estimator. Standard error forumlae for the population size 
estimates can be found in Stoklosa et al. (2011).
}

\value{ Output gives population size estimates and standard errors. }

\references{Stoklosa, Hwang, Wu and Huggins. (2011). Heterogeneous 
capture-recapture models with covariates: A partial likelihood approach for 
closed populations. \emph{Biometrics}.}

\author{Jakub Stoklosa, Wenhan Hwang and Richard Huggins.}

\keyword{models}