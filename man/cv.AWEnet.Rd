\name{cv.AWEnet}
\alias{cv.AWEnet}
\title{Computes K-fold cross validated error curve for AEnet and WEnet}
\description{This function computes the K-fold cross validation estimates.}
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
  \item{K}{number of folds.}
  \item{fraction}{abscissa values at which CV curve should be computed. This is the fraction of the saturated |beta|. The default value is seq(from = 0, to = 1, length =100).}
  \item{plot.it}{if T then plot will be showed. Default is T.}
  \item{se}{include standard error bands.}
  \item{AEnet}{if T then the results are based on adaptive elastic net otherwise based on weighted elastic net.}
   \item{all.folds}{null.}
}
\details{This function computes the K-fold cross validation, cross validation error, cross validation mean squared error.}

\value{An "index" object is returned with a CV curve. The index includes
\item{lambda2}{as AEnetCC.aft}
\item{cv}{the CV curve at each value of index}
\item{cv.mse}{the mean square error of the CV curve}
\item{cv.error}{the standard error of the CV curve}

}

\references{
Khan and Shaw (2015) imputeYn: Imputing the last largest censored
observation/observations under weighted least squares. R package version
1.3, \url{https://cran.r-project.org/package=imputeYn}.

Khan and Shaw (2015). Variable Selection for Survival Data with a Class of Adaptive Elastic Net Techniques. Statistics and Computing (published online; DOI: 10.1007/s11222-015-9555-8). Also available in \url{http://arxiv.org/abs/1312.2079}.
}
\author{Hasinur Rahaman Khan and Ewart Shaw}

\seealso{cv.AWEnetCC}
\examples{
#For full data typically used for AFT models (using imputeYn (2015) package).
dat<-data(n=100, p=10, r=0, b1=c(rep(5,5),rep(0,5)), sig=1, Cper=0)

#This needs to run for generating weights of the observations
l<-mrbj(cbind(dat$y, dat$delta) ~ dat$x, mcsize=100, trace=FALSE, gehanonly=FALSE)

#cv.AWEnet: Cross validation of Adaptive elastic net
wt<-l$enet
\dontrun{cv1 <-cv.AWEnet(dat$x, dat$y, dat$delta, weight=wt, lambda2=0.001, maxit=10,
plot.it = T, AEnet=T)}
\dontrun{cv1$index[which.min(cv1$cv)]}

#cv.AWEnet: Cross validation of weighted elastic net
\dontrun{l<-mrbj(cbind(dat$y, dat$delta) ~ dat$x, mcsize=100, trace=FALSE, gehanonly=TRUE)}
\dontrun{wt<-l$gehansd}
\dontrun{cv2 <-cv.AWEnet(dat$x, dat$y, dat$delta, weight=wt, lambda2=0.001, 
maxit=10, plot.it = T, AEnet=F)}
\dontrun{cv2$index[which.min(cv2$cv)]}
}
\keyword{CV}
\keyword{AEnet}
\keyword{WEnet}
