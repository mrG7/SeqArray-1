\name{seqClose}
\alias{seqClose}
\title{Close a Sequencing GDS File}
\description{
	Close a Sequencing GDS file
}
\usage{
seqClose(gdsfile)
}
\arguments{
	\item{gdsfile}{a \code{\link{SeqVarGDSClass}} object}
}
\value{
	Return an object of class \code{\link[gdsfmt]{gds.class}}.
}
\details{
	It is strongly suggested to call \code{seqClose} instead of \code{\link[gdsfmt]{closefn.gds}}.
}

\author{Xiuwen Zheng}
\seealso{
	\code{\link{seqOpen}}
}

\examples{
gds.fn <- seqExampleFileName("gds")
# or gds.fn <- "C:/YourFolder/Your_GDS_File.gds"

# open the GDS file
gdsfile <- seqOpen(gds.fn)

# display the contents of the GDS file in a hierarchical structure
gdsfile

# close the GDS file
seqClose(gdsfile)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
