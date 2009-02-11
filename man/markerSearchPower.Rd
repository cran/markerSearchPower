\name{markerSearchPower}
\alias{markerSearchPower}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function of power calculation }
\description{
  Analytically calculate statistical power for three model selection methods in genome-wide association studies: marginal selection, exhaustive search, and forward selection. Power definitions are (A) the probability of identifying exactly the true model (in marginal selection, it is the probability of detecting both true markers); and (b) the probability of detecting at least one of the true markers. 
}
\usage{
markerSearchPower(mainEff1, mainEff2, epistasisEff, noiseSD, alleleFreq1, alleleFreq2, 
                 alleleFreq3 = 0.5, DetectionN, obsN, markerN, 
                 strategy = c("marginal", "exhaustive", "forward"), 
                 powerDef = c("both", "either"), samplePointN = 20000,
		 decompMethod = c("svd", "chol", "eigen"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mainEff1}{ The main effect of associated marker 1. It is represented by b1 in the true genetic model. See Details below.}
  \item{mainEff2}{ The main effect of associated marker 2. It is represented by b2 in the true genetic model. See Details below.}
  \item{epistasisEff}{ The epistatic effect (interaction) between markers 1 and 2. It is represented by b3 in the true genetic model. }
  \item{noiseSD}{ The standard deviation of noise term. }
  \item{alleleFreq1}{ The frequency of the effective allele at marker 1.}
  \item{alleleFreq2}{ The frequency of the effective allele at marker 2.}
  \item{alleleFreq3}{ Allele frequency of other non-associated markers. This value has little effect on power calculation.}
  \item{DetectionN}{ The number of top models selected. For power of detecting both associated markers in marginal selection, \code{DetectionN >= 2}.}
  \item{obsN}{ Sample size. }
  \item{markerN}{ Number of markers genptyped. }
  \item{strategy}{ The stratgy of marker searching: "marginal", "exhaustive", or "forward". Default is "marginal". }
  \item{powerDef}{ Definition of power. Either "all" indicates the probability of identifying exactly the true model (in marginal selection, it is the probability of detecting both true markers); or "either" indicates the probability of detecting at least one of the true markers. }
  \item{samplePointN}{ Number of random points for Monte Carlo integration. Default is 20,000. }
  \item{decompMethod}{ Matrix decomposition used to determine the matrix root of sigma, possible methods are singular value decomposition ("svd", default), eigenvalue decomposition ("eigen"), and Cholesky decomposition ("chol"). This argument is only useful for forward selection identifying either true SNP. }
}
\details{
  The power calculation is based on a hypothetical true underlying genetic model of a two biallelic marker system for a quantitative trait: y = b0 + b1*X1 + b2*X2 + b3*X1*X2 + noise. Noise has normal distribution with mean zero and standard deviation \code{noiseSD}. X1 and X2 are the genotypes of markers 1 and 2, respectively. Genotype values are based on the assumption of additive allelic effect: it is equal to the number of copies of the effective allele. For the main covariates, the signs of coefficients b1 and b2 illustrate the direction of effective alleles of markers 1 and 2. The coefficient b3 for the interaction term represents the direction and magnitude of epistasis. The power calculation for binary trait is under development. 
}
\value{
 The power of a given marker search strategy under a given definition either (A) or (B).
}
\references{Zheyang Wu and Hongyu Zhao (2009)
  \emph{Statistical Power of Model Selection Strategies for Genome-Wide Association Studies}. Submitted.
 }
\author{ 
Zheyang Wu and Hongyu Zhao 
(Maintainer: Zheyang Wu <zheyang.wu@yale.edu>)
}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

markerSearchPower (mainEff1=0.1, mainEff2=0.1, epistasisEff=2.4, noiseSD=sqrt(3), 
     alleleFreq1=0.7, alleleFreq2=0.7, DetectionN=2, obsN=100, markerN=300);
markerSearchPower (mainEff1=0.1, mainEff2=0.1, epistasisEff=2.4, noiseSD=sqrt(3), 
     alleleFreq1=0.7, alleleFreq2=0.7, DetectionN=1, obsN=100, markerN=300, 
     powerDef = "either");
markerSearchPower (mainEff1=0.1, mainEff2=0.1, epistasisEff=1.4, noiseSD=sqrt(3), 
     alleleFreq1=0.7, alleleFreq2=0.7, DetectionN=1, obsN=100, markerN=300, 
     strategy = "exhaustive", powerDef = "either");
markerSearchPower (mainEff1=0.1, mainEff2=0.1, epistasisEff=1.4, noiseSD=sqrt(3), 
     alleleFreq1=0.7, alleleFreq2=0.7, DetectionN=1, obsN=100, markerN=300, 
     strategy = "exhaustive", powerDef = "both");
markerSearchPower (mainEff1=0.1, mainEff2=0.1, epistasisEff=1.4, noiseSD=sqrt(3), 
     alleleFreq1=0.7, alleleFreq2=0.7, DetectionN=5, obsN=100, markerN=300, 
     strategy = "forward", powerDef = "both");
markerSearchPower (mainEff1=0.1, mainEff2=0.1, epistasisEff=1.4, noiseSD=sqrt(3), 
     alleleFreq1=0.7, alleleFreq2=0.7, DetectionN=5, obsN=100, markerN=300, 
     strategy = "forward", powerDef = "either");

%## The function is currently defined as
%function (mainEff1, mainEff2, epistasisEff, noiseSD, alleleFreq1, 
%    alleleFreq2, alleleFreq3 = 0.5, DetectionN, obsN, markerN, 
%    strategy = c("marginal", "exhaustive", "forward"), powerDef = c("both", 
%        "either"), samplePointN = 20000, decompMethod = c("svd", 
%        "chol", "eigen")) 
%{
%    b1 = mainEff1
%    b2 = mainEff2
%    b3 = epistasisEff
%    sigm = noiseSD
%    p1 = alleleFreq1
%    p2 = alleleFreq2
%    p3 = alleleFreq3
%    R = DetectionN
%    n = obsN
%    p = markerN
%    meanT1 = meanT1Func(b1, b2, b3, p1, p2, sigm, n)
%    meanT2 = meanT1Func(b2, b1, b3, p2, p1, sigm, n)
%    varT1 = varT1Func(b1, b2, b3, p1, p2, sigm)
%    varT2 = varT1Func(b2, b1, b3, p2, p1, sigm)
%    covT1T2 = covT1T2Func(b1, b2, b3, p1, p2, sigm)
%    covMT1T2 = matrix(c(varT1, covT1T2, covT1T2, varT2), ncol = 2)
%    meanT12 = meanT12Func(b1, b2, b3, p1, p2, sigm, n)
%    varT12 = varT12Func(b1, b2, b3, p1, p2, sigm)
%    covT12T1 = covT12T1Func(b1, b2, b3, p1, p2, sigm)
%    covT12T2 = covT12T1Func(b2, b1, b3, p2, p1, sigm)
%    covMT12s = matrix(c(varT12, covT12T1, covT12T2, covT12T1, 
%        varT1, covT1T2, covT12T2, covT1T2, varT2), ncol = 3)
%    meanF3l1 = meanF3l1Func(b1, b2, b3, p1, p2, p3, sigm)
%    varF3l1 = varF3l1Func(b1, b2, b3, p1, p2, p3, sigm)
%    c3l1 = varF3l1/(2 * meanF3l1)
%    d3l1 = 2 * meanF3l1^2/varF3l1
%    meanF3l2 = meanF3l1Func(b2, b1, b3, p2, p1, p3, sigm)
%    varF3l2 = varF3l1Func(b2, b1, b3, p2, p1, p3, sigm)
%    c3l2 = varF3l2/(2 * meanF3l2)
%    d3l2 = 2 * meanF3l2^2/varF3l2
%    meanT1l3 = meanT1l3Func(b1, b2, b3, p1, p2, sigm, n)
%    meanT2l3 = meanT1l3Func(b2, b1, b3, p2, p1, sigm, n)
%    varT1l3 = varT1l3Func(b1, b2, b3, p1, p2, sigm)
%    varT2l3 = varT1l3Func(b2, b1, b3, p2, p1, sigm)
%    covT1T1l3 = covT1T1l3Func(b1, b2, b3, p1, p2, sigm)
%    covT2T2l3 = covT1T1l3Func(b2, b1, b3, p2, p1, sigm)
%    covT1T2l3 = covT1T2l3Func(b1, b2, b3, p1, p2, sigm)
%    covT2T1l3 = covT1T2l3Func(b2, b1, b3, p2, p1, sigm)
%    covT1l3T2l3 = covT1l3T2l3Func(b1, b2, b3, p1, p2, sigm)
%    covMTvector = matrix(c(varT1, covT1T2, covT1T1l3, covT1T2l3, 
%        covT1T2, varT2, covT2T1l3, covT2T2l3, covT1T1l3, covT2T1l3, 
%        varT1l3, covT1l3T2l3, covT1T2l3, covT2T2l3, covT1l3T2l3, 
%        varT2l3), ncol = 4)
%    if (strategy == "marginal" && powerDef == "both") {
%        R = DetectionN - 1
%        return(power.marginal(meanT1, meanT2, covMT1T2, R, p, 
%            TRUE, samplePointN))
%    }
%    if (strategy == "marginal" && powerDef == "either") 
%        return(power.marginal(meanT1, meanT2, covMT1T2, R, p, 
%            FALSE, samplePointN))
%    if (strategy == "exhaustive" && powerDef == "both") 
%        return(power.exhaustive.both(meanT12, meanT1, meanT2, 
%            covMT12s, c3l1, d3l1, c3l2, d3l2, R, p, n, samplePointN))
%    if (strategy == "exhaustive" && powerDef == "either") 
%        return(power.exhaustive.atleast1(meanT12, meanT1, meanT2, 
%            covMT12s, c3l1, d3l1, c3l2, d3l2, R, p, n, samplePointN))
%    if (strategy == "forward" && powerDef == "both") 
%        return(power.forward.both(meanT12, meanT1, meanT2, covMT12s, 
%            c3l1, d3l1, c3l2, d3l2, R, p, n, samplePointN))
%    if (strategy == "forward" && powerDef == "either") 
%        return(power.forward.atleast1(meanT1, meanT2, meanT1l3, 
%            meanT2l3, covMTvector, R, p, samplePointN, decompMethod[1]))
%    cat(paste("Power calculation for finding ", powerDef, " true markers with "%, 
%        strategy, " search strategy failed", sep = ""))
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
