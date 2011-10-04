\name{PL.popN}

\alias{PL.popN}

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

\usage{PL.popN(model, tau, D)}

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

\value{Output gives population size estimates and standard errors.}

\references{Stoklosa, Hwang, Wu and Huggins. (2011). Heterogeneous 
capture-recapture models with covariates: A partial likelihood approach for 
closed populations. \emph{Biometrics}.}

\author{Jakub Stoklosa, Wenhan Hwang and Richard Huggins.}

\examples{
## Several modelling examples using the Mountain Pygmy Possum capture-recapture 
## data set as given in PL.popN.

library(PL.popN)

data(possum) ## Load data.

tau <- 5     ## No. of capture occasions.
possum$y     ## No. of times captured.
possum$t1     ## Time of first capture.
possum$x     ## Covariate (body weight).
D.poss <- nrow(possum)   ## No. of individuals captured at least once.

h <- tau - possum$t1     ## No. of trials.
R <- (possum$y - 1)/h   ## Observed proportions.
R[is.na(R)] <- 0
possum <- transform(possum, h = h, R = R)

## Quadratic GLM.
poss.glm <- glm(R ~ x + I(x^2), family = binomial,
                weights = h, data = possum)
PL.popN(poss.glm, tau, D.poss)

## GAM.
library(mgcv)
poss.gam <- gam(R ~ s(x, bs = "cr"), family = binomial, 
                weights = h, data = possum, gamma = 1.4)
PL.popN(poss.gam, tau, D.poss)
  
## Plots.
pred <- predict.gam(poss.gam, type = "response", se.fit = TRUE)
plot.mat <- cbind(pred$fit, pred$fit - 2*pred$se.fit,
                  pred$fit + 2*pred$se.fit)
par(las = 1, font.lab = 2)
matplot(sort(possum$x), plot.mat[order(possum$x),], 
        col = c(1,2,2), lty = c(1,2,2), xlab = "Body weight",
        ylab = "Fitted capture probabilities", type = 'l',
        font = 2, cex.axis = 1.5, lwd = 3, cex.lab = 1.8)

## Simex.
library(simex)
x.poss <- possum$x
poss.glm <- glm(R ~ x.poss, family = binomial,
                x = TRUE, y = TRUE, weights = h)
poss.simex <- simex(poss.glm, measurement.error = 0.5, 
                    SIMEXvariable = "x.poss")
PL.popN(poss.simex, tau, D.poss)

## GLMM.
library(MASS)
nID <- 1 : D.poss
poss.glmm <- glmmPQL(R ~ x.poss, random = ~ 1 | nID, 
                     family = binomial, verbose = FALSE)
PL.popN(poss.glmm, tau, D.poss)

## Local polynomial model.
library(sm)
y.poss <- possum$y-1   ## Observed counts.
last.cap <- which(possum$t1==tau)
y.poss <- y.poss[-last.cap]
x.poss <- x.poss[-last.cap]
N.poss <- (tau-possum$t1)[-last.cap]

poss.lp <- sm.binomial(scale(x.poss), y.poss, N.poss,
                       h = 1, nbins = 20, eval.points = sort(scale(possum$x)))
PL.popN(poss.lp, tau, D.poss)
} 

\keyword{models}