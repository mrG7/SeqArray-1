SeqVarGDSClass <- function(gdsobj, ...) {
    new("SeqVarGDSClass", gdsobj, ...)
}

setMethod("granges",
          "SeqVarGDSClass",
          function(x, ...) {
            variant.id <- seqGetData(x, "variant.id")
            chromosome <- seqGetData(x, "chromosome")
            position <- seqGetData(x, "position")
            gr <- GRanges(seqnames=chromosome,
                          ranges=IRanges(start=position, end=position),
                          ...)
            names(gr) <- variant.id
            gr
          })

setMethod("ref",
          "SeqVarGDSClass",
          function(x) {
            a <- seqGetData(x, "allele")
            ref <- substr(a, 1, regexpr(",", a)-1)
            DNAStringSet(ref)
          })

setMethod("alt",
          "SeqVarGDSClass",
          function(x) {
            a <- seqGetData(x, "allele")
            alt <- substr(a, regexpr(",", a)+1, nchar(a))
            do.call(DNAStringSetList, strsplit(alt, ",", fixed=TRUE))
          })

