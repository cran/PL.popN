\name{PL.popN0}

\alias{PL.popN0}

\title{Closed population size estimation from capture-recapture data using a 
partial likelihood with covariates.}

\description{\code{PL.popN0} gives closed population size estimates and the 
standard errors for capture-recature partial likelihood models when capture 
propabilities depended on covariates. Models recognized by \code{PL.popN0} must 
come in the form of the following logistic regression models: Generalized Linear
Models (GLM), Generalized Additive Models (GAM), Generalized Linear Mixed Model
(GLMM) or Simulation Extrapolation (SIMEX) models. The number of capture 
occasions must also be specified.
}

\usage{PL.popN0(model, tau)}

\arguments{
\item{model}{ an object of model type, i.e. model objects from either 
\code{glm, gam, glmmPQL, simex}. }
\item{tau}{ the number of capture occasions in the experiment. }
}

\details{The structure of the partial likelihood allows one to construct 
logistic regression models and analyse the data parametrically by using 
\code{glm()}, non-parametric (semi-parametric) smoothing by using \code{gam()}, 
mixed models by using \code{glmmPQL} and measurement error by usng \code{simex}. 
Population size estimates are obtained by extracting fitted capture 
probabilities from the models, which are then employed in a Horwitz-Thompson 
estimator. Standard error forumlae for the population size estimates can be 
found in Stoklosa et al. (To appear).
}

\value{ Output gives population size estimates and standard errors. }

\references{Stoklosa, Hwang, Wu and Huggins. (In Press). Heterogeneous 
capture-recapture models with covariates: A partial likelihood approach for 
closed populations. Biometrics.}

\author{Jakub Stoklosa, Wenhan Hwang and Richard Huggins.}

\keyword{models}
