setOldClass("gds.class")
setClass("SeqVarGDSClass", contains="gds.class")

setValidity("SeqVarGDSClass",
    function(object) {
        n <- index.gdsn(object, "description", silent=TRUE)
        if (is.null(n))
            return("Description variable must exist!")
        tmp <- get.attr.gdsn(n)
        if (!("sequence.variant.format" %in% names(tmp)))
            return("Not a sequencing-variant GDS file!")

        var.names <- ls.gdsn(object)
        if (!all(c("sample.id", "variant.id", "position",
                "chromosome", "allele", "genotype") %in% var.names))
            return("sample.id, variant.id, position, chromosome, allele, and genotype are required variables.")
        TRUE
    })
