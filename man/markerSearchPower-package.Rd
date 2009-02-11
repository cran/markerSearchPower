\name{markerSearchPower-package}
\alias{markerSearchPower-package}
%\alias{markerSearchPower}
\docType{package}
\title{
Package of power calculation for marker detection strategies 
}
\description{
This package calculates statistical power of detecting associated markers based on one of the model selection strategies: marginal selection, exhaustive search, or forward selection. With assumed genetic effects (including interaction effect), allele frequencies, noise level, sample size, number of genotyped markers, and control level (i.e. number of markers or models intended to select), this package provides fast and accurate consultation on power of different model selection strategies. It helps researchers to decide a more efficient way for marker detection in genome-wide association studies.  
}
\details{
\tabular{ll}{
Package: \tab markerSearchPower\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-02-11\cr
License: \tab GPL-2\cr
%LazyLoad: \tab yes\cr
}


This package requires three other R packages: mvtnorm, adapt, corpcor.


Installation from local files: 
For windows:
Download the .zip file into your local hard drive.
Open Rgui (R graphical interface), go to menu "Packages or install packages
from local (zip) files ...", then select the provided zip file.  Go to menu 
"Packages or load package..." to load the package into R for use.  
For Linux or Unix:
Download the .tar.gz file into your local hard drive.
Install the R package from the source file using the following commands:
first go to the directory where the downloaded file is stored and then use
"R CMD INSTALL <filename>.tar.gz". For more options please see R help on
"INSTALL".


Very easy to use. Only one external function. Choose proper parameters and methods for function \code{markerSearchPower}, which returns a power value.
}
\author{
Zheyang Wu and Hongyu Zhao

Maintainer: Zheyang Wu <zheyang.wu@yale.edu>

}
\references{
Zheyang Wu and Hongyu Zhao (2009)
  \emph{Statistical Power of Model Selection Strategies for Genome-Wide Association Studies}. Submitted.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
markerSearchPower (mainEff1=0.1, mainEff2=0.1, epistasisEff=2.4, noiseSD=sqrt(3), alleleFreq1=0.7, 
            alleleFreq2=0.7, strategy = "marginal", powerDef = "both", DetectionN=10,
            obsN=1000, markerN=300000);

}
