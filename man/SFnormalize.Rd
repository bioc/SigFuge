\name{SFnormalize}
\alias{SFnormalize}
\title{SigFuge normalize read counts}
\usage{
  SFnormalize(data, flag = 1)
}
\arguments{
  \item{data}{a \eqn{d \times n}{d x n} matrix of read
  counts at \eqn{d} positions for \eqn{n} samples.}

  \item{flag}{a \eqn{n \times 1}{n x 1} logical vector of
  samples flagged as low expression. If \code{flag == 1},
  default low expression cutoffs are used. If \code{flag ==
  0}, no samples are flagged as low expression (equivalent
  to setting flag = zeros(n, 1)).}
}
\value{
  \code{SFnormalize} returns a list containing: \itemize{
  \item data.norm a \eqn{d \times (n-m)}{d x (n-m)} matrix
  of normalized read counts where \eqn{m} is the number 
  of low expression samples.  \item flag a \eqn{n \times
  1}{n x 1} logical vector of flagged samples.}
}
\description{
  Function for normalizing read count data as specified in
  the SigFuge method. The normalization procedure is
  applied prior to SigFuge clustering to remove the effect
  of sample-locus specific expression from the analysis.
  This allows the method to identify clusters based on
  expression patterns across the genomic locus. It is
  recommended to flag and remove low expression samples
  from the normalization and analysis since their shapes may
  be overwhelmed by noise. A threshold based method for
  identifying low expression samples is included in the
  function, but users may also specify their own flags
  for low expression samples.
}
\examples{
data(geneDepth)
depthnorm <- SFnormalize(geneDepth, flag = 1)
}
\author{
  Patrick Kimes <pkimes@live.unc.edu>
}

