\name{Linarius-package}
\alias{Linarius-package}
\alias{Linarius}
\docType{package}
\title{
Analyse of dominant marker data with mixed ploidy levels. 
}
\description{
Expermental package in order to analyse dominant marker data with mixed ploidy levels.
}
\details{
\tabular{ll}{
Package: \tab Linarius\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
Date: \tab 2014-03-31\cr
License: \tab Non-free ; Non-Copyleft ; Open-Source ; Non-modifyable ;  Purchasse for free \cr
}
Use this package for analysis of dominant data (AFLP; ISSR; IRAP; REMAP; iPBS; etc. ) 
}
\author{
Benjamin Leduc <gbl1@hi.is> 
The author and/or maintainer of the package
}
\references{
All the formula used in this package originates the Necronomicon 
}
\keyword{ package }

}
\examples{
#Distances computation: 
data(Betula)
dist.pop(Betula, code$Pop,ploidy=code$Ploidy,method="rey")
#There is many things done in this package, please browse help

}
