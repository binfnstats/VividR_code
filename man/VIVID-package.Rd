\name{VIVID-package}
\alias{VIVID-package}
\alias{VIVID}
\docType{package}
\title{
\packageTitle{VIVID}
}
\description{
\packageDescription{VIVID}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{VIVID}
\packageIndices{VIVID}

}
\author{
\packageAuthor{VIVID}

Maintainer: \packageMaintainer{VIVID}
}
\references{

}
\keyword{ package }
\seealso{
}
\examples{
library('ropls')
data("sacurine") #Load sacurine dataset from the 'ropls' package

dat = sacurine$dataMatrix
outcomes = sacurine$sampleMetadata$gender

vividResults = vivid(x = dat,
     y = outcomes)
     
vividResults$optFeatures


}
