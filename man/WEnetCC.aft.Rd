\name{WEnetCC.aft}
\alias{WEnetCC.aft}
\title{Weighted elastic net with censoring constraints for censored data based on AFT models}
\description{This function is used to fit AFT models using weighted elastic net with censorig constraints (Khan and Shaw, 2013). This is an extension of the weighted elastic approach that allows the censoring constraints to be implemented into the optimization framework. The method works for both cases of datasets high-dimensional where the number of variables (p) is greater than the number of subjects (n) and low-dimensional.}

\usage{WEnetCC.aft(X, Y, delta, weight, C, s, lambda2)}

\arguments{
  \item{X}{covariate matrix under study, particularly for AFT modelling. The order of matrix covariate is typically n by p.}
  \item{Y}{typically the logarithmic of the survival time under AFT models. Otherwise survival time.}
  \item{delta}{status. it includes value 1 for uncensored and value 0 for censored subject.}
  \item{weight}{vector of observation weights.}
  \item{C}{this is a positive value that accounts for the penalties of violations of constraints. C is typically allowed to take values in a grid such as (0, 0.5, 1, 1.5, ..., 10).}
  \item{s}{this is the optimal equivalent specification for lasso in terms of fraction of the L1 norm. This is obtained from the AEnet.aft function.}
  \item{lambda2}{regularization parameter for the L2 norm of the coefficients. This is typically assumed to take values in a relatively small grid, say (0, 0.5, 1.0, 1.5, 2.0, ..., 5).}
}

\details{This method use the same optimal pair of (lambda1, lambda2) as found in AEnet.aft. Then C is typically allowed to take values in a grid such as (0, 0.5, 1.0, 1.5, 2.0,...,5), and the optimal value for C obtained by 5-fold cross-validation. Here C typically depends upon how stringently one wants the model to satisfy the censoring constraints compared to how good is the prediction for uncensored data.}

\value{An object of type numeric is returned that provides the coefficient estimates.}

\references{
Hong, D. and Zhang, F. (2010). Weighted Elastic Net Model for Mass
Spectrometry Imaging Processing. Mathematical Modelling of Natural
Phenomena, 5, 115-133.

Khan and Shaw (2013). Variable Selection for Survival Data with A Class of Adaptive Elastic Net Techniques. CRiSM working paper, Department of Statistics, University of Warwick, UK.
}

\author{Hasinur Rahaman Khan and Ewart Shaw}

\seealso{cv.AWEnetCC}

\examples{
#For full data typically used for AFT models (using imputeYn package)
dat<-data(n=100, p=10, r=0, b1=c(rep(5,5),rep(0,5)), sig=1, Cper=0)

#This needs to run for generating weights of the observations
l<-mrbj(cbind(dat$y, dat$delta) ~ dat$x, mcsize=10, trace=FALSE, gehanonly=TRUE)

#WEnetCC.aft: weighted elastic net with censoring constraints
wt<-round(l$gehansd)
ft<-WEnetCC.aft(dat$x, dat$y, dat$delta, weight=wt, C=1, s=959596, lambda2=0.5)
\donttest{ft}
}

\keyword{WEnetCC.aft}
\keyword{AFT models}
\keyword{lars}
