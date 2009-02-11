\name{markerSearchPower-internal}
\title{Internal functions}
\alias{.required}
\description{
Internal markerSearchPower functions. 
}
\usage{
CDF.F3l1(x, c, d)
CDF.F4l3(x)
CDF.F34(x)
CDF.T3(x)
covT1l3T2l3Func(b1, b2, b3, p1, p2, sigm)
covT1T1l3Func(b1, b2, b3, p1, p2, sigm)
covT1T2Func(b1, b2, b3, p1, p2, sigm)
covT1T2l3Func(b1, b2, b3, p1, p2, sigm)
covT12T1Func(b1, b2, b3, p1, p2, sigm)
meanF3l1Func(b1, b2, b3, p1, p2, p3, sigm)
meanT1Func(b1, b2, b3, p1, p2, sigm, n)
meanT1l3Func(b1, b2, b3, p1, p2, sigm, n)
meanT12Func(b1, b2, b3, p1, p2, sigm, n)
PDF.F34(x)
power.exhaustive.atleast1(meanT12, meanT1, meanT2, covMT12s, c3l1, 
                          d3l1, c3l2, d3l2, R, p, n, samplePointN)
power.exhaustive.both(meanT12, meanT1, meanT2, covMT12s, c3l1, d3l1, 
                     c3l2, d3l2, R, p, n, samplePointN)
power.forward.atleast1(meanT1, meanT2, meanT1l3, meanT2l3, covMTvector, 
                       R, p, samplePointN, decompMethod)
power.forward.both(meanT12, meanT1, meanT2, covMT12s, c3l1, d3l1, 
                   c3l2, d3l2, R, p, n, samplePointN)
power.marginal(meanT1, meanT2, covMT1T2, R, p, isBoth, samplePointN)
tripleSet(r, p, F34N)
tripleSet2(r, p1, p2, F34N)
varF3l1Func(b1, b2, b3, p1, p2, p3, sigm)
varT1Func(b1, b2, b3, p1, p2, sigm)
varT1l3Func(b1, b2, b3, p1, p2, sigm)
varT12Func(b1, b2, b3, p1, p2, sigm)
}
%\arguments{
%  \item{x}{ vector of quantiles. }
%  \item{c, d}{ parameters for weighted chi-sqare distribution. }
%}

\details{
  These functions are not for use at user level. 
%\item{CDF.T3}is for cumulative distribution function of T-test statistics, which are asymptotically standard normal.
}
\value{
  \code{CDF.F3l1, CDF.F4l3, CDF.F34, CDF.T3} give discribution functions.
  \code{PDF.F34} gives density functions.
  \code{meanF3l1Func, meanT1Func, meanT1l3Func, meanT12Func} give means for relevant asymptotic distributions.  
  \code{varF3l1Func, varT1Func, varT1l3Func, varT12Func} give variances for relevant asymptotic distributions.  
  \code{covT1l3T2l3Func, covT1T1l3Func, covT1T2Func, covT1T2l3Func, covT12T1Func} give covariances for relevant asymptotic distributions. 
  \code{power.xx} give calculated power values.  
  \code{tripleSet, tripleSet2} return all possible distributions in three groups for finite samples. 
}
\references{ 
Zheyang Wu and Hongyu Zhao (2009)
  \emph{Statistical Power of Model Selection Strategies for Genome-Wide Association Studies}. Submitted 
}
\author{ Zheyang Wu }

\keyword{internal}
