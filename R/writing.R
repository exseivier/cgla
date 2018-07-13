#############
# LIBRARIES #
#############

# // NO LIBRARIES NEEDED.
# // LIBRARIES ARE LOADED IN extract.R

###########
# OBJECTS #
###########

# // NO OBJECTS

###########
# HELPERS #
###########

# /*
#  * Writes a bounch of sequences stored in DNAStringSet
#  * allocated in a list object.
#  * This function requires a DNAStringSetList: a list of
#  * DNAStringSet objects, and a character vector with the
#  * name of the output path.
#  * The names of the output files will be taken from the
#  * DNAStringSetList names.
#  */
write.DNAStringSetList <- function(path=NULL, dssl) {
	if (is.null(path)) {
	    path <- "."
	}
	for (name in names(dssl)){
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

###########
# METHODS #
###########

# // NO METHODS
