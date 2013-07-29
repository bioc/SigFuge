\name{SFpval}
\alias{SFpval}
\title{Calculate SigFuge \eqn{p}-value}
\usage{
  SFpval(data, normalize = 1, flag = 1)
}
\arguments{
  \item{data}{a \eqn{d \times n}{d x n} matrix of read
  counts at \eqn{d} positions for \eqn{n} samples.}

  \item{normalize}{a \eqn{n \times 1}{n x 1} logical vector
  of flagged samples.}

  \item{flag}{a \eqn{n \times 1}{n x 1} logical vector of
  samples flagged as low expression. If \code{flag == 1},
  default low expression cutoffs are applied to
  \code{data}. If \code{flag == 0}, no samples are flagged
  as low expression (equivalent to setting \code{flag =
  zeros(n,1)}).}
}
\value{
  \code{SFpval} returns an object of class
  \code{\link{sigclust-class}}. Avaliable slots are
  described in detail in the \code{\link{sigclust}}
  package. Primarily, we make use of \code{@pvalnorm}.
}
\description{
  Function for computing significance of clustering
  \eqn{p}-value. \eqn{p}-value is obtained from
  \code{\link{sigclust}}, a simulation based procedure for
  testing significance of clustering in high dimension low
  sample size (HDLSS) data.

  The SigClust hypothesis test is given: \itemize{ \item
  H0: data generated from single Gaussian \item H1: data
  not generated from single Gaussian }
}
\examples{
data(geneDepth)
SFout <- SFpval(geneDepth, normalize = 1, flag = 1)
SFout@pvalnorm
}
\author{
  Patrick Kimes <pkimes@live.unc.edu>
}

