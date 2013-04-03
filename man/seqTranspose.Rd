\name{seqTranspose}
\alias{seqTranspose}
\title{Transpose Data Array}
\description{
	Transpose data array or matrix for possibly higher-speed access
}
\usage{
seqTranspose(gdsfile, var.name, compress=NULL, verbose=TRUE)
}
\arguments{
	\item{gdsfile}{a \code{\link[gdsfmt]{gds.class}} object in the \link{gdsfmt} package}
	\item{var.name}{the variable name with '/' as a separator}
	\item{compress}{the compression option used in \code{\link[gdsfmt]{add.gdsn}}; or
		determine automatically if \code{NULL}}
	\item{verbose}{if TRUE, show information}
}
\value{
	None.
}
\details{
	It is designed for possibly higher-speed access. More details will be provided
in the future version.
}

\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link[SeqArray]{seqGetData}}, \code{\link[SeqArray]{seqApply}}
}

\examples{
# the file of VCF
(vcf.fn <- seqExampleFileName("vcf"))
# or vcf.fn <- "C:/YourFolder/Your_VCF_File.vcf"

# convert
seqVCF2GDS(vcf.fn, "tmp.gds")

# list the structure of GDS variables
f <- seqOpen("tmp.gds", FALSE)
f

seqTranspose(f, "genotype/data")
f

# the original array
index.gdsn(f, "genotype/data")
# the transposed array
index.gdsn(f, "genotype/~data")

# close
seqClose(f)
unlink("tmp.gds")
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}