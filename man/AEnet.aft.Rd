\name{AEnet.aft}
\alias{AEnet.aft}
\title{Adaptive elastic net for censored data based on AFT models}
\description{This function is used to fit AFT models using adaptive elastic net.}
\usage{AEnet.aft(X, Y, delta, weight, lambda2, maxit = 10)}
\arguments{
  \item{X}{covariate matrix under study, particularly for AFT modelling. The order of matrix covariate is typically n by p.}
  \item{Y}{typically the logarithmic of the survival time under AFT models. Otherwise survival time.}
  \item{delta}{status. it includes value 1 for uncensored and value 0 for censored subject.}
  \item{weight}{vector of observation weights. Weight is based on initial estimator that is obtained from elastic net on the weighted data (see Enet.wls function).}
  \item{lambda2}{regularization parameter for the L2 norm of the coefficients. This is typically assumed to take values in a relatively small grid, say (0, 0.5, 1.0, 1.5, 2.0, ..., 5).}
  \item{maxit}{an optional bound for the number of steps to be taken. Default is 10.}
}
\details{
This function is used to fit AFT models using adaptive elastic net approach (Khan and Shaw, 2015). The method works for both cases of datasets high-dimensional where the number of variables (p) is greater than the number of subjects (n) and low-dimensional.

The adaptive elastic net is transformed into an adaptive lasso type problem in an augmented space and then is solved using the lars algorithm. This is an extension of the work Ghosh (2007) for censored data. To find the optimal value for the tuning parameters lambda11 and lambda2, first lambda2 is typically assumed to take values in a relatively small grid, say (0, 0.5, 1.0, 1.5, 2.0, ..., 5). For each lambda2, the lars algorithm produces the entire solution path. This gives the optimal equivalent specification for lasso in terms of fraction of the L1 norm.}
\value{A "AEnet.aft" object is returned. It includes
\item{beta}{coefficient estimates of the covariates}
\item{mu}{means of Y}
\item{meanx}{means of the columns of X}
\item{normx}{normalized value of the columns of X}
\item{type}{lasso}
}

\references{
Ghosh, S. (2007). Adaptive Elastic Net: An Improvement of Elastic Net
to achieve Oracle Properties. Technical Reports, Indiana University-
Purdue University, Indianapolis, USA. PR no. 07-01.

Khan and Shaw (2015) imputeYn: Imputing the last largest censored
observation/observations under weighted least squares. R package version
1.3, \url{https://cran.r-project.org/package=imputeYn}.

Khan and Shaw (2015). Variable Selection for Survival Data with a Class of Adaptive Elastic Net Techniques. Statistics and Computing (published online; DOI: 10.1007/s11222-015-9555-8). Also available in \url{http://arxiv.org/abs/1312.2079}.

}
\author{Hasinur Rahaman Khan and Ewart Shaw}

\seealso{cv.AWEnet}
\examples{
#For full data typically used for AFT models (using imputeYn (2015) package).
dat<-data(n=100, p=10, r=0, b1=c(rep(5,5),rep(0,5)), sig=1, Cper=0)

#This needs to run for generating weights of the observations
l<-mrbj(cbind(dat$y, dat$delta) ~ dat$x, mcsize=100, trace=FALSE, gehanonly=FALSE)

#AEnet.aft: adaptive elastic net
wt<-round(l$enet)
ft.1<-AEnet.aft(dat$x, dat$y, dat$delta, weight=wt, lambda2=1, maxit=10)
\donttest{ft.1}
}
\keyword{Adaptive elastic net}
\keyword{AFT models}
\keyword{lars}
