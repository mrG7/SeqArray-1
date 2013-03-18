\name{seqMerge}
\alias{seqMerge}
\title{Merge Multiple Sequencing GDS Files}
\description{
	Merge multiple sequencing GDS files
}
\usage{
seqMerge(gds.fn, out.fn, compress.option = seqCompress.Option(),
	verbose = TRUE)
}
\arguments{
	\item{gds.fn}{the file names of multiple GDS files}
	\item{out.fn}{the output file name}
	\item{compress.option}{specify the compression options, by default
		\code{\link{seqCompress.Option}}}
	\item{verbose}{if TRUE, show information}
}
\value{
	None.
}
\details{
	The current implementation of \code{seqMerge} extracts and merges the genotypic
data only without any annotation. Users can specify multiple VCF files in
\code{\link[SeqArray]{seqVCF2GDS}} to export a single GDS file.
}

\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link[SeqArray]{seqVCF2GDS}}
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
