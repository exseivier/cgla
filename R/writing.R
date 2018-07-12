# LIBRARIES

# OBJECTS


# HELPERS

write.DNAStringSetList <- function(path=NULL, dssl) {
	for (name in names(dssl)){
		if (is.null(path)) {
			path <- "."
		}
		output <- strsplit(name, "/")[[1]]
		output <- output[length(output)]
		output <- strsplit(output, "\\.")[[1]]
		output <- output[1:(length(output) - 1)]
		output <- paste(output, collapse=".")
		output <- paste(output, "_interSeqs.fsa", sep="")
		output <- paste(path, output, sep="/")
		writeXStringSet(dssl[[name]], output, format="fasta")
	}
}
# METHODS
