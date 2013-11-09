\name{AdapEnetClass-package}
\docType{package}
\title{A class of adaptive elastic net methods for censored data}
\description{The package provides new approaches to variable selection for censored data, based on AFT models optimized using regularized weighted least squares. Namely, a weighted elastic net, an adaptive elastic net, and two of their extensions by adding censoring observations as constraints into their model optimization frameworks.
}
\details{
The accelerated failure time (AFT) models have proved useful in many contexts, though heavy censoring (as for example in cancer survival) and high dimensionality (as for example in microarray data) cause difficulties for model fitting and model selection. The package provide new approaches to variable selection for censored data, based on AFT models optimized using regularized weighted least squares. The regularized technique uses a mixture of L1 and L2 norm penalties under two elastic net type approaches. The approaches extend the original approaches proposed by Ghosh (2007), and Hong and Zhang (2010). This package also provides two extended approaches by adding censoring observations as constraints into their model optimization frameworks.

\tabular{ll}{
Package: \tab AdapEnetClass\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-11-09\cr
License: \tab GPL-2\cr
Depends: \tab imputeYn, glmnet, lars\cr}
}
\author{
Hasinur Rahaman Khan and Ewart Shaw
Maintainer: Hasinur Rahaman Khan <hasinurkhan@gmail.com>
}

\references{
Ghosh, S. (2007). Adaptive Elastic Net: An Improvement of Elastic Net
to achieve Oracle Properties. Technical Reports, Indiana University-
Purdue University, Indianapolis, USA. PR no. 07-01.

Hong, D. and Zhang, F. (2010). Weighted Elastic Net Model for Mass
Spectrometry Imaging Processing. Mathematical Modelling of Natural
Phenomena, 5, 115-133.

Jin, Z., Lin, D. Y. and Ying, Z. (2006). On least-squares regression with censored data. Biometrika, 93, 147-161.

Khan and Shaw (2013). Variable Selection for Survival Data with A Class of Adaptive Elastic Net Techniques. CRiSM working paper, Department of Statistics, University of Warwick, UK.

Khan and Shaw (2013b). Variable Selection with The Modified Buckley-James Method and The Dantzig Selector for High-dimensional Survival Data. Conference Proceeding, 59th World Statistics Congress, Hong Kong, China
}
\keyword{lars}
\keyword{imputeYn}

\examples{
#For full data typically used for AFT models (using imputeYn package)
dat<-data(n=100, p=10, r=0, b1=c(rep(5,5),rep(0,5)), sig=1, Cper=0)

#mrbj: modified resampling based buckley-james
l<-mrbj(cbind(dat$y, dat$delta) ~ dat$x, mcsize=100, trace=FALSE, gehanonly=FALSE)

#AEnet.aft: adaptive elastic net
wt<-round(l$enet)
ft.1<-AEnet.aft(dat$x, dat$y, dat$delta, weight=wt, lambda2=1, maxit=10)

#AEnetCC.aft: adaptive elastic net with censoring constraints
\donttest{ft.1cc<-AEnetCC.aft(dat$x, dat$y, dat$delta, weight=wt, C=1, 
s = 0.959596, lambda2=0.5)}

#WEnet.aft: weighted elastic net
#mrbj: modified resampling based buckley-james
\donttest{l<-mrbj(cbind(dat$y, dat$delta) ~ dat$x, mcsize=100, trace=FALSE, gehanonly=TRUE)}
\donttest{wt<-l$gehansd}
\donttest{ft.2<-WEnet.aft(dat$x, dat$y, dat$delta, weight=wt, lambda2=0.5, maxit=10)}

#WEnetCC.aft: weighted elastic net with censoring constraints
\donttest{ft.2cc<-WEnetCC.aft(dat$x, dat$y, dat$delta, weight=wt, C=1, s = 1, lambda2=0.5)}
}
