\name{plotObsEst}
\alias{plotObsEst}
\title{Pairwise scatter plots of the survival times}
\description{This function generates pairwise scatter plots of the observed and predicted survival times.}
\usage{plotObsEst(yObs, yEst, delta, xlab = NULL, ylab = NULL, title = NULL, 
legendplot = TRUE, legendpos = "topleft", maxvalue = NULL, minvalue = NULL)}
\arguments{
  \item{yObs}{observed survival times.}
  \item{yEst}{predicted survival times.}
  \item{delta}{status. it includes value 1 for uncensored and value 0 for censored subject.}
  \item{xlab}{x-axis title. Default is NULL.}
  \item{ylab}{y-axis title. Default is NULL.}
  \item{title}{main title of the plot.}
  \item{legendplot}{if TRUE, plot the legend.}
  \item{legendpos}{position of the legend in the plot.}
  \item{maxvalue}{maximum value of the yObs and yEst. Default is Null.}
  \item{minvalue}{minimum value of the yObs and yEst. Default is Null.}
}
\value{This provides a scatter plot of the survival times.}
\references{
Khan and Shaw (2013). Variable Selection for Survival Data with a Class of Adaptive Elastic Net Techniques. CRiSM working paper, Department of Statistics, University of Warwick, UK, No. 13-17. Also available in \url{http://arxiv.org/abs/1312.2079}.}

\examples{
#For a hypothetical data
y<-c(12,33,22,34,11,23)
delta<-c(1,0,0,1,0,1)
yEst<-c(11,30,24,30,6,13)

#plotObsEst: scatter plotting of the pairwise survival times 
\donttest{plot<-plotObsEst(y, yEst, delta, xlab = NULL, ylab = NULL, 
title = Predicted versus Observed Survival times, legendplot = TRUE,
legendpos = "topleft", maxvalue = NULL, minvalue = NULL)}
}

