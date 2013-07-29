\name{SFfigure}
\alias{SFfigure}
\title{Plot expression as curves}
\usage{
  SFfigure(data, locusname, annot = c(), flip.fig = 1,
    label.exon = 1, print.n = 1, data.labels = 0,
    label.colors = c(), flag = 1, lplots = 2,
    savestr = c(), titlestr = c(), pval = 1)
}
\arguments{
  \item{data}{a \eqn{d \times n}{d x n} matrix or
  data.frame of read counts at \eqn{d} base positions for
  \eqn{n} samples.}

  \item{locusname}{a character string specifying gene or
  locus name to be used in figure title.}

  \item{annot}{a \code{GRanges} object or \code{data.frame}
  including annotation information for locus, including:
  \itemize{ \item \code{start} start of contiguous genomic
  regions \item \code{end} end of contiguous genomic
  regions \item \code{seqname} chromosome name for genomic
  region \item \code{strand} strandedness of sequence }}

  \item{flip.fig}{an indicator whether to flip the plotting
  direction of the locus if \code{strand == "-"} when
  annotation information is provided.}

  \item{label.exon}{an indicator whether to print the exon
  boundaries to the figure.}

  \item{print.n}{an indicator whether to print cluster
  sizes.}

  \item{data.labels}{a \eqn{n \times 1}{n x 1} vector of
  class labels to use instead of calcuating SigFuge labels}

  \item{label.colors}{a \eqn{K \times 3}{K x 3} matrix of
  RBG colors specifying cluster colors for \eqn{K} clusters.
  \code{ggplot2} default colors are used if not specified. 
  If using SigFuge default labels, \eqn{K=3} even if no 
  low expression samples are flagged.}

  \item{flag}{a \eqn{n \times 1}{n x 1} logical vector of
  samples flagged as low expression. If \code{flag == 1},
  default low expression cutoffs are used. If \code{flag ==
  0}, no samples are flagged as low expression (equivalent
  to setting flag = zeros(n, 1)).}

  \item{lplots}{a specification of which figures to output
  \itemize{ \item \code{1}: curves in single
  panel, random colors  \item \code{2}: curves in single panel, 
  colored by cluster \item \code{3}: curves in \eqn{K} panels, 
  separated and colored by cluster  \item \code{4}: curves in 
  \eqn{n} panels, colored by cluster (single sample per panel)
  \item \code{5}: cluster medians in single panel, colored by 
  cluster }}

  \item{savestr}{a string specifying the file name for
  resulting figures. Extensions can also be specified in
  \code{savestr}. If no extension is specified figures will
  be saved as pdfs. If \code{length(lplots) > 1}, figures
  will be saved as \code{paste0(savestr,"_x")} for \code{x}
  in \code{lplots} with the appropirate extension. If no
  \code{savestr} is specified, function will return a list
  containing the created \code{ggplot} objects.}

  \item{titlestr}{a string specifying figure title. If
  unspecified, default is 
  \code{titlestr=}\code{paste(locusname," locus, SigFuge
  analysis")}.}

  \item{pval}{an indicator whether the \code{SFpval}
  should be computed. If \code{pval == 1}, the
  p-value is added to the title, i.e.
  (\code{titlestr=paste0(titlestr, ", p-value = ", p)}).}
}
\value{
  \code{SFfigure} returns a figure that is saved to the
  current working directory if a \code{savestr} is
  specified. Else, a list containing the plots is returned.
}
\description{
  Function for producing various figures corresponding to
  the SigFuge functional data approach to studying RNA-seq
  data as expression curves along base positions. The primary 
  input for the function is a read count matrix and 
  \code{GRanges}. The default behavior is to identify 
  clusters based on applying \code{SFlabels} to a normalized
  version of the data produced by \code{SFnormalize}. 
  If specified, the function will compute a p-value for the 
  significance of the labels by calling the \code{SFpval}
  function.
}
\examples{
# load data
data(geneAnnot)
data(geneDepth)

# only use first 50 samples
mdata <- geneDepth[,1:50]

# make plot
locusname <- "CDKN2A"
SFfigure(mdata, locusname, geneAnnot, flag=1,
 lplots=3, savestr=paste0(locusname,".pdf"), titlestr="CDKN2A locus, LUSC samples",
 pval=1)

mySFs <- SFfigure(mdata, locusname, geneAnnot, flag=1,
            lplots=1, savestr=c(), titlestr="CDKN2A locus, LUSC samples not saved",
            pval=0)
mySFs$plot1
}
\author{
  Patrick Kimes <pkimes@live.unc.edu>
}

