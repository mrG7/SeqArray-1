\name{seqVCF.Header}
\alias{seqVCF.Header}
\title{Parse the header of a VCF file}
\description{
	Parse the header of a VCF file
}
\usage{
seqVCF.Header(vcf.fn)
}
\arguments{
	\item{vcf.fn}{the file name of VCF}
}
\value{
	Return a list (with class name "seqvcf.header.class"):
	\item{fileformat}{the file format}
	\item{info}{the ID description in the INFO field}
	\item{filter}{the ID description in the FILTER field}
	\item{format}{the ID description in the FORMAT field}
	\item{alt}{the ID description in the ALT field}
	\item{contig}{the description in the contig field}
	\item{assembly}{the link of assembly}
	\item{header}{the other header lines}
	\item{num.ploidy}{the number of ploidy, two for humans}
}
\details{
	The ID description contains four columns:
	ID -- variable name;
	Number -- the number of elements, see the webpage of the 1000 Genomes Project;
	Type -- data type;
	Description -- a variable description.
}

\references{
	\url{http://www.1000genomes.org/wiki/Analysis/Variant Call Format/vcf-variant-call-format-version-41}

	The variant call format and VCFtools.
	Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE,
	Lunter G, Marth GT, Sherry ST, McVean G, Durbin R; 1000 Genomes Project Analysis Group.
	Bioinformatics. 2011 Aug 1;27(15):2156-8. Epub 2011 Jun 7.
}
\author{Xiuwen Zheng \email{zhengx@u.washington.edu}}
\seealso{
	\code{\link[SeqArray]{seqVCF.SampID}}, \code{\link[SeqArray]{seqVCF2GDS}}
}

\examples{
(vcf.fn <- seqExampleFileName("vcf"))
# or vcf.fn <- "C:/YourFolder/Your_VCF_File.vcf"

# get sample id
seqVCF.Header(vcf.fn)
}

\keyword{gds}
\keyword{sequencing}
\keyword{genetics}