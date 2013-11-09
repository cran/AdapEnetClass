\name{mrbj}
\alias{mrbj}
\title{Modified resampling based Buckley-James method}
\description{This function is based on an alternative method to the Buckley-James method based on the accelerated failure time models as discussed in Jin et al.(2006) that provides an approximating approach for the consistent root of the estimating equation since the general Buckley-James method does not guarantee the convergence of its iterative process (Khan and Shaw, 2013b).}
\usage{
mrbj(formula, data, subset, trace = FALSE, gehanonly = FALSE, cov = FALSE, 
na.action = na.exclude, residue = FALSE, mcsize = 100)
}
\arguments{
  \item{formula}{specifies a model to be fitted. The response and covariates of the model are separated by a ~ operator. The response, on the left side of ~, should be a
Surv object with two columns, of which the first column is the survival time or
censored time and the second column is the censoring indicator. The covariates
or predictors X, on the right side of ~, should be columns with the same length
as Surv object.}
  \item{data}{a data frame which contains the Surv objects and covariates.}
  \item{subset}{specifies subset of the original data frame that should be used for the model fit.}
  \item{trace}{takes logical values TRUE or FALSE. If it is set to be TRUE, then the summary of every iteration will be kept. The default is FALSE.}
  \item{gehanonly}{takes logical values T or F. If gehanonly=T, only Gehan estimator will be calculated and otherwise, estimator will be calculated using elastic net applied to weighted data. The default is gehanonly=F.}
  \item{cov}{takes logical values T or F. If cov=T, the covariance matrices of the Gehen estimator will be printed. The default is cov=F.}
  \item{na.action}{takes values na.exclude or na.fail. The default is na.exclude, which deletes the observations with missing values. The other choice is na.fail, which returns an
error if any missing values are found.}
 \item{residue}{Default is FALSE.}
 \item{mcsize}{specifies the resampling number. The default is 500.}
}

\details{Modified resampling based Buckley James method is an extension of the work proposed by Jin et al. (2006). Here, if n (observations)>p (number of covariates under AFT models) then Gehan estimator is used an an initial estimator in the Buckley-James iterative process to produce consistent estimator. For n<p, the elastic net on the weighted data is used and bootsrap is used to generate standard error of the estimates.}
\value{The Gehan estimator, the standard error of the Gehan estimator, the Z score and the p-value for testing the hypothesis of beta=0 based on Gehan estimation. The elastic net estimator, the standard error of the elastic net estimator, The covariance matrices of the Gehan estimator, when cov is set to be T.}
\references{
Jin, Z., Lin, D. Y. and Ying, Z. (2006). On least-squares regression with censored data. Biometrika, 93, 147-161.

Khan and Shaw (2013b). Variable Selection with The Modified Buckley-James Method and The Dantzig Selector for High-dimensional Survival Data. Conference Proceeding, 59th World Statistics Congress, Hong Kong, China

}
\author{Hasinur Rahaman Khan and Ewart Shaw}


\seealso{lss package}

\examples{
#For full data typically used for AFT models (using imputeYn package)
dat<-data(n=100, p=10, r=0, b1=c(rep(5,5),rep(0,5)), sig=1, Cper=0)

#mrbj: modified resampling based Buckley-James method. 
#When gehanonly=T, it give both the Gehan and the elastic net  estimates 
\donttest{fit1<-mrbj(cbind(dat$y, dat$delta) ~ dat$x, mcsize=100, trace=FALSE, gehanonly=TRUE)}
\donttest{fit1}

#mrbj: modified resampling based Buckley-James method. 
#Only for Gehan estimates
fit2<-mrbj(cbind(dat$y, dat$delta) ~ dat$x, mcsize=100, trace=FALSE, gehanonly=FALSE)
\donttest{fit2}
}

\keyword{Buckley-James method}
\keyword{Elastic net}
