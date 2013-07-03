\name{seqVCF2GDS}
\alias{seqVCF2GDS}
\title{Reformat VCF files}
\description{
	Reformat Variant Call Format (VCF) files
}
\usage{
seqVCF2GDS(vcf.fn, out.fn, header = NULL, genotype.var.name = "GT",
	compress.option = seqCompress.Option(),
	cvt.raise.err = TRUE, verbose = TRUE)
}
\arguments{
	\item{vcf.fn}{the file name(s) of VCF format}
	\item{out.fn}{the file name of output GDS file}
	\item{header}{if NULL, \code{header} is set to be \code{\link{seqVCF.Header}(vcf.fn)}}
	\item{genotype.var.name}{the ID for genotypic data in the FORMAT column;
		"GT" by default, VCFv4.0}
	\item{compress.option}{specify the compression options, by default
		\code{\link{seqCompress.Option}}}
	\item{cvt.raise.err}{\code{TRUE}: throw an error if numeric conversion fails;
		\code{FALSE}: get missing value if numeric conversion fails}
	\item{verbose}{if TRUE, show information}
}
\value{
	Return the file name with an absolute path.
}
\details{
	GDS -- Genomic Data Structures used for storing genetic array-oriented data,
		and the file format used in the \link{gdsfmt} package.

	VCF -- The Variant Call Format (VCF), which is a generic format for storing DNA
		polymorphism data such as SNPs, insertions, deletions and structural variants,
		together with rich annotations.

	If there are more than one files in \code{vcf.fn}, \code{seqVCF2GDS} will
merge all dataset together once they all contain the same samples. It is useful to
combine genetic data if VCF data are divided by chromosomes.
}
\references{
	The variant call format and VCFtools.
	Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE,
	Lunter G, Marth GT, Sherry ST, McVean G, Durbin R; 1000 Genomes Project Analysis Group.
	Bioinformatics. 2011 Aug 1;27(15):2156-8. Epub 2011 Jun 7.

	\url{http://corearray.sourceforge.net/}
}

\author{Xiuwen Zheng}
\seealso{
	\code{\link{seqVCF.Header}}, \code{\link{seqCompress.Option}},
	\code{\link{seqGDS2VCF}}
}

\examples{
# the file of VCF
vcf.fn <- seqExampleFileName("vcf")
# or vcf.fn <- "C:/YourFolder/Your_VCF_File.vcf"

# convert
seqVCF2GDS(vcf.fn, "tmp.gds")

# display
(f <- seqOpen("tmp.gds"))

# close the GDS file
seqClose(f)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}
