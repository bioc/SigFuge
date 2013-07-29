\name{SFlabels}
\alias{SFlabels}
\title{Calculate SigFuge labels}
\usage{
  SFlabels(normData)
}
\arguments{
  \item{normData}{a list containing \itemize{ \item
  \code{data.norm} a \eqn{d \times (n-m)}{d x (n-m)} matrix
  of normalized read counts at \eqn{d} positions for
  \eqn{(n-m)} samples where \eqn{n} is the total number of 
  samples and \eqn{m} is the number of low expression 
  samples. \item \code{flag} a \eqn{n \times 1}{n x
  1} logical vector of flagged samples with \eqn{\sum
  \texttt{flag} = m}{sum(flag) = m}. }}
}
\value{
  \code{SFlabels} returns a \eqn{n \times 1}{n x 1}
  vector of class labels.
}
\description{
  Function for producing vector of SigFuge labels using 2-means 
  clustering on non-low expression normalized data and combining
  with low expression flags. Typically, \code{\link{SFlabels}} 
  is used by passing output from \code{\link{SFnormalize}}.
}
\examples{
data(geneDepth)
normalizedData <- SFnormalize(geneDepth)
labels <- SFlabels(normalizedData)
}
\author{
  Patrick Kimes <pkimes@live.unc.edu>
}

