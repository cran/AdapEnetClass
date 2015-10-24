\name{Enet.wls}
\alias{Enet.wls}
\title{Elastic net based on the weighted data}
\description{This function provides the estimates of the initial estimators.}
\arguments{
  \item{X}{covariate matrix under study, particularly for AFT modelling. The order of matrix covariate is typically n by p.}
  \item{Y}{typically the logarithmic of the survival time under AFT models. Otherwise survival time.}
  \item{delta}{status. it includes value 1 for uncensored and value 0 for censored subject.}
}
\details{This function provides the estimates of the initial estimators that are used to calculate weights for observations under adaptive and weighted elastic net approaches. Here simple elastic net method is applied to the data weighted by the Kaplan-Meier weights (Khan and Shaw, 2015). The elastic net estimates are obtained using glmnet.}

\value{
\item{beta}{coefficient estimates of the covariates}
\item{fit}{an object that gives df (The number of nonzero coefficients for each value of lambda), Dev (deviance) abd Lambda (The actual sequence of lambda values used)}
}

\references{
Khan and Shaw (2015) imputeYn: Imputing the last largest censored
observation/observations under weighted least squares. R package version
1.3, \url{https://cran.r-project.org/package=imputeYn}.

Khan and Shaw (2015). Variable Selection for Survival Data with a Class of Adaptive Elastic Net Techniques. Statistics and Computing (published online; DOI: 10.1007/s11222-015-9555-8). Also available in \url{http://arxiv.org/abs/1312.2079}.
}
\author{Hasinur Rahaman Khan and Ewart Shaw}

\examples{
#For full data typically used for AFT models (using imputeYn (2015) package).
dat<-data(n=100, p=10, r=0, b1=c(rep(5,5),rep(0,5)), sig=1, Cper=0)

#Enet.wls: Elastic net on the weighted data
enet<-Enet.wls(dat$x, dat$y, dat$delta)
\donttest{enet}
}

\keyword{elastic net}
\keyword{Kaplan-Meier weights}
