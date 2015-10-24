\name{cv.AWEnetCC}
\alias{cv.AWEnetCC}
\title{Computes K-fold cross validated error curve for AEnetCC and WEnetCC}
\description{This function computes the K-fold cross validation estimates.}
\usage{
cv.AWEnetCC(X, Y, delta, weight, kFold = 10, C, s, lambda2, AEnetCC = T)
}
\arguments{
  \item{X}{covariate matrix under study, particularly for AFT modelling. The order of matrix covariate is typically n by p.}
  \item{Y}{typically the logarithmic of the survival time under AFT models. Otherwise survival time.}
  \item{delta}{status. it includes value 1 for uncensored and value 0 for censored subject.}
  \item{weight}{vector of observation weights. Weight is based on initial estimator that is obtained from elastic net on the weighted data (see Enet.wls function) or from Gehan estimator (see mrbj function).}
  \item{kFold}{number of folds.}
  \item{C}{this is a positive value that accounts for the penalties of violations of constraints. C is typically allowed to take values in a grid such as (0, 0.5, 1, 1.5, ..., 10).}
  \item{s}{this is the optimal equivalent specification for lasso in terms of fraction of the L1 norm. This is obtained from the AEnet.aft function}.
  \item{lambda2}{regularization parameter for the L2 norm of the coefficients. This is typically assumed to take values in a relatively small grid.}
  \item{AEnetCC}{If T then the results are based on adaptive elastic net with censoring constraints otherwise based on the weighted elastic net with censoring constraints.}
}

\details{The function gives the K-fold cross validation, cross validation error, cross validation mean squared error.}
\value{
  \item{beta}{shows coefficient estimates of the covariates.}
  \item{betavar}{variance of the coefficient estimates.}
  \item{cvscore}{a CV score based on the CV error. This is basically the sum of squared residuals of uncensored data multiplied by the Kaplan-Meier weights (Khan and Shaw, 2015).} 
}
\references{
Khan and Shaw (2015) imputeYn: Imputing the last largest censored
observation/observations under weighted least squares. R package version
1.3, \url{https://cran.r-project.org/package=imputeYn}.

Khan and Shaw (2015). Variable Selection for Survival Data with a Class of Adaptive Elastic Net Techniques. Statistics and Computing (published online; DOI: 10.1007/s11222-015-9555-8). Also available in \url{http://arxiv.org/abs/1312.2079}.
}
\seealso{cv.AWEnet}

\examples{
#For full data typically used for AFT models (using imputeYn (2015) package)
dat<-data(n=100, p=10, r=0, b1=c(rep(5,5),rep(0,5)), sig=1, Cper=0)

#This needs to run for generating weights of the observations
l<-mrbj(cbind(dat$y, dat$delta) ~ dat$x, mcsize=100, trace=FALSE, gehanonly=FALSE)

#cv.AWEnetCC: Cross validation of Adaptive elastic net with censoring constraints
wt<-l$enet
cv1cc<-cv.AWEnetCC(dat$x, dat$y, dat$delta, weight=wt, kFold = 10, C=1.2, s=0.88,
 lambda2=0.001, AEnetCC=TRUE)

#cv.AWEnetCC: Cross validation of weighted elastic net with censoring constraints
\dontrun{l<-mrbj(cbind(dat$y, dat$delta) ~ dat$x, mcsize=100, trace=FALSE, gehanonly=TRUE)}
\dontrun{wt<-l$gehansd}
\dontrun{cv1cc<-cv.AWEnetCC(dat$x, dat$y, dat$delta, weight=wt, kFold = 10, C=1.2, s=0.88,
 lambda2=0.001, AEnetCC=F)}
}
\keyword{CV}
\keyword{AEnetCC}
\keyword{WEnetCC}
