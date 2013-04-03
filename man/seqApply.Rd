\name{seqApply}
\alias{seqApply}
\title{Apply Functions Over Array Margins}
\description{
	Returns a vector or list of values obtained by applying a function to margins
of arrays or matrices
}
\usage{
seqApply(gdsfile, var.name, FUN,
	margin = c("by.variant"),
	as.is = c("list", "integer", "double", "character", "none"),
	var.index = FALSE, ...)
}
\arguments{
	\item{gdsfile}{a \code{\link[gdsfmt]{gds.class}} object in the \link{gdsfmt} package}
	\item{var.name}{the variable name(s), see details}
	\item{FUN}{the function to be applied}
	\item{margin}{giving the dimension which the function will be applied
		over. E.g., for a matrix 1 indicates rows, 2 indicates columns}
	\item{as.is}{returned value: a list, an integer vector, etc}
	\item{var.index}{if TRUE, add an argument to the user-defined function
		\code{FUN} like \code{FUN(index, x, ...)} where \code{index} is an index of
		variant starting from 1 if \code{margin = "by.variant"};
		otherwise call \code{FUN(x, ...)}}
	\item{...}{optional arguments to \code{FUN}}
}
\details{
	The variable name should be \code{"sample.id"}, \code{"variant.id"},
\code{"position"}, \code{"chromosome"}, \code{"allele"}, \code{"annotation/id"},
\code{"annotation/qual"}, \code{"annotation/filter"},
\code{"annotation/info/VARIABLE_NAME"}, or \code{"annotation/format/VARIABLE_NAME"}.

	The algorithm is highly optimized by blocking the computations to exploit the
high-speed memory instead of disk.
}
\value{
	A vector or list of values.
}

\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link[SeqArray]{seqSetFilter}},
	\code{\link[SeqArray]{seqGetData}},
	\code{\link[SeqArray]{seqParallel}}
}

\examples{
# the file of GDS
gds.fn <- seqExampleFileName("gds")
# or gds.fn <- "C:/YourFolder/Your_GDS_File.gds"

# display
(f <- seqOpen(gds.fn))

# get 'sample.id
(samp.id <- seqGetData(f, "sample.id"))
# "NA06984" "NA06985" "NA06986" ...

# get 'variant.id'
head(variant.id <- seqGetData(f, "variant.id"))


# set sample and variant filters
set.seed(100)
seqSetFilter(f, sample.id=samp.id[c(2,4,6,8,10)],
	variant.id=sample(variant.id, 10))

# read multiple variables variant by variant
seqApply(f, c(geno="genotype", phase="phase", qual="annotation/id"),
	FUN=function(x) print(x), as.is="none")

# get the numbers of alleles per variant
seqApply(f, "allele",
	FUN=function(x) length(unlist(strsplit(x,","))), as.is="integer")


################################################################
# with an index of variant

seqApply(f, c(geno="genotype", phase="phase", qual="annotation/id"),
	FUN=function(index, x) { print(index); print(x); index },
	as.is="integer", var.index=TRUE)
# it is as the same as
which(seqGetFilter(f)$variant.sel)



################################################################
# reset sample and variant filters
seqSetFilter(f)

# calculate the frequency of reference allele,
#   a faster version could be obtained by C coding
af <- seqApply(f, "genotype", FUN=function(x) mean(x==0, na.rm=TRUE), as.is="double")
length(af)
summary(af)


# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}