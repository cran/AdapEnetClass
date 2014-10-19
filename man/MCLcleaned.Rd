\name{MCLcleaned}
\docType{data}
\alias{MCLcleaned}
\title{Mantle cell lymphoma cleaned data}
\description{This is a cleaned gene expression Mantle cell lymphoma dataset that contains expression values of 574 cDNA elements as used in Khan and Shaw (2013). The full dataset is available
at http://llmpp.nih.gov/MCL/ as used in Rosenwald et al. (2003). Total 92 patients were classified as having MCL, based on established morphologic and immunophenotypic criteria. Survival times of 64 patients were
available and the remaining 28 patients were censored. The data do not provide any further relevant covariates for MCL patients. Use as "data(MCLcleaned)".}
\usage{data(MCLcleaned)}

\format{ It contains total 577 columns namely, ID: array identifier for 92 patients, time: time of follow-up in year, cens: pateint status at follow up (1=death, 0=censored), and the remaining 574 columns are the gene expression value of 574 cDNA denoted by X37,..,X8893 which are the measurements of the log2 CY5/CY3 expression ratio for each feature. The dummy gene id X.... can be matched with the original gene UNIQID as stated at \url{http://llmpp.nih.gov/MCL/}.
}


\source{The full dataset is available at \url{http://llmpp.nih.gov/MCL/}.}


\references{
Khan and Shaw (2013). Variable Selection for Survival Data with a Class of Adaptive Elastic Net Techniques. CRiSM working paper, Department of Statistics, University of Warwick, UK, No. 13-17. Also available in \url{http://arxiv.org/abs/1312.2079}.

Rosenwald A., Wright G., Wiestner A., Chan W., Connors J., Campo E., Gascoyne
R., Grogan T., Muller Hermelink H., Smeland E., Chiorazzi M., Giltnane J., Hurt E., Zhao H., Averett
L., Henrickson S., Yang L., Powell J., Wilson W., Jaffe E., Simon R., Klausner R., Montserrat E., Bosch
F., Greiner T., Weisenburger D., Sanger W., Dave B., Lynch J., Vose J., Armitage J., Fisher R., Miller T.,
LeBlanc M., Ott G., Kvaloy S., Holte H., Delabie J., Staudt L. (2003): The proliferation gene expression signature
is a quantitative integrator of oncogenic events that predicts survival in mantle cell lymphoma. Cancer Cell,
3, 185-197.
}
\author{Hasinur Rahaman Khan and Ewart Shaw}

\examples{
#Enet.wls: Elastic net on the weighted data Using MCLcleaned data attached with this package
\donttest{data(MCLcleaned)}
\donttest{GEval<-MCLcleaned[4:577]}
\donttest{enet<-Enet.wls(GEval, log(MCLcleaned@time), MCLcleaned@cens)}
\donttest{enet}
}

\keyword{dataset}

