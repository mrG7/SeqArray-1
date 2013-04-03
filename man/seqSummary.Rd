\name{seqSummary}
\alias{seqSummary}
\title{Get the summary of a GDS file}
\description{
	Get the summary of a sequencing GDS file
}
\usage{
seqSummary(gdsfile, varname=NULL, check=c("check", "full.check", "none"), verbose=TRUE)
}
\arguments{
	\item{gdsfile}{a file name, or a \code{\link[gdsfmt]{gds.class}} object in
		the \link{gdsfmt} package}
	\item{varname}{if \code{NULL}, check the whole GDS file; or a character specifying
		variable name, and return a description of that variable. See details.}
	\item{check}{should be one of "check", "full.check", "none"}
	\item{verbose}{if \code{TRUE}, display information}
}
\value{
	If \code{varname = NULL}, then return a list:
		\item{filename}{the file name}
		\item{sequence.variant.format}{the sequencing format in GDS}
		\item{num.of.sample}{the number of samples}
		\item{num.of.variant}{the number of variants}
		\item{info}{the description of the INFO field: \code{var.name}, \code{number},
			\code{type} and \code{description}}
		\item{format}{the description of the FORMAT field: \code{var.name}, \code{number},
			\code{type} and \code{description}}

	If \code{varname = "genotype" or "phase"}, then return a list:
		\item{dim}{dim[1] -- ploidy, dim[2] -- the number of samples,
			dim[3] -- the number of variants}
		\item{seldim}{seldim[1] -- the number of selected samples,
			seldim[2] -- the number of selected variants}
}
\details{
	If \code{check = "check"}, this function performs regular checking: dimensions
of variables, etc. If \code{check = "full.check"}, it performs more checking: unique
sample id, unique variant id, whether genotypic data are in a valid range or not, etc.
}

\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link[SeqArray]{seqGetData}}, \code{\link[SeqArray]{seqApply}}
}

\examples{
(gds.fn <- seqExampleFileName("gds"))
# or gds.fn <- "C:/YourFolder/Your_GDS_File.gds"

seqSummary(gds.fn)

seqSummary(gds.fn, "genotype")


######################################################
# display
f <- seqOpen(gds.fn)

# get 'sample.id
samp.id <- seqGetData(f, "sample.id")
# get 'variant.id'
variant.id <- seqGetData(f, "variant.id")

# set sample and variant filters
seqSetFilter(f, sample.id=samp.id[c(2,4,6,8,10)])
set.seed(100)
seqSetFilter(f, variant.id=sample(variant.id, 10))

seqSummary(f, "genotype")

# close a GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}