\name{cv.AWEnet}
\alias{cv.AWEnet}
\title{Computes K-fold cross validated error curve for AEnet and WEnet}
\description{This function computes the K-fold cross validation, cross validation error, cross validation mean squared error.}
\usage{
cv.AWEnet(X, Y, delta, weight, lambda2, maxit, K = 10, fraction = seq(from = 0, 
to = 1, length = 100), plot.it = F, se = TRUE, AEnet = T, all.folds = NULL)}
\arguments{
  \item{X}{covariate matrix under study, particularly for AFT modelling. The order of matrix covariate is typically n by p.}
  \item{Y}{typically the logarithmic of the survival time under AFT models. Otherwise survival time.}
  \item{delta}{status. it includes value 1 for uncensored and value 0 for censored subject.}
  \item{weight}{vector of observation weights. Weight is based on initial estimator that is obtained from elastic net on the weighted data (see Enet.wls function) or from Gehan estimator (see mrbj function).}
  \item{lambda2}{regularization parameter for the L2 norm of the coefficients. This is typically assumed to take values in a relatively small grid.}
  \item{maxit}{an optional bound for the number of steps to be taken. Default is 10.}
  \item{K}{Number of folds.}
  \item{fraction}{Abscissa values at which CV curve should be computed. This is the fraction of the saturated |beta|. The default value is seq(from = 0, to = 1, length =100).}
  \item{plot.it}{If T then plot will be showed. Default is T.}
  \item{se}{Include standard error bands.}
  \item{AEnet}{If T then the results are based on adaptive elastic net otherwise based on weighted elastic net.}
   \item{all.folds}{NULL.}
}

\value{An "index" object is returned with a CV curve. The index includes
\item{lambda2}{as above}
\item{cv}{The CV curve at each value of index}
\item{cv.mse}{The mean square error of the CV curve}
\item{cv.error}{The standard error of the CV curve}

}

\references{
Khan and Shaw (2013). Variable Selection for Survival Data with A Class of Adaptive Elastic Net Techniques. CRiSM working paper, Department of Statistics, University of Warwick, UK.
}
\author{Hasinur Rahaman Khan and Ewart Shaw}

\seealso{cv.AWEnetCC}
\examples{
#For full data typically used for AFT models (using imputeYn package).
dat<-data(n=100, p=10, r=0, b1=c(rep(5,5),rep(0,5)), sig=1, Cper=0)

#This needs to run for generating weights of the observations
l<-mrbj(cbind(dat$y, dat$delta) ~ dat$x, mcsize=100, trace=FALSE, gehanonly=FALSE)

#cv.AWEnet: Cross validation of Adaptive elastic net
wt<-l$enet
\donttest{cv1 <-cv.AWEnet(dat$x, dat$y, dat$delta, weight=wt, lambda2=0.001, maxit=10,
plot.it = T, AEnet=T)}
\donttest{cv1$index[which.min(cv1$cv)]}

#cv.AWEnet: Cross validation of weighted elastic net
\donttest{l<-mrbj(cbind(dat$y, dat$delta) ~ dat$x, mcsize=100, trace=FALSE, gehanonly=TRUE)}
\donttest{wt<-l$gehansd}
\donttest{cv2 <-cv.AWEnet(dat$x, dat$y, dat$delta, weight=wt, lambda2=0.001, 
maxit=10, plot.it = T, AEnet=F)}
\donttest{cv2$index[which.min(cv2$cv)]}
}
\keyword{CV}
\keyword{AEnet}
\keyword{WEnet}
