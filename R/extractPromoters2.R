#############
# LIBRARIES #
#############

library(S4Vectors)
library(Biostrings)

###########
# OBJECTS #
###########

# /* S4 object to store data from tab annotation file */
setClass("pocket", representation(chrname="character", name="character", start="numeric", end="numeric", strand="character"))

# /* S4 objet used to store the sequence and its coordinates, and name and chromosome name */
setClass("pocketOfSeq", representation(chrname="character", name="character", seq="DNAStringSet", coords="numeric", strand="character"))

# /* S4 object used as list of pocket or pocketOfSeq objects */
setClass("listOfPockets", representation(list="list"), validity = function(object){
	if(length(object@list) > 3){
		stop("listOfPockets is higher than 3")
	}
})

####################
# OBJECTS WRAPPERS #
####################

NEW.POCKET <- function(chrname, name, start, end, strand){
	obj <- new("pocket", chrname=chrname, name=name, start=start, end=end, strand=strand)
	obj
}

###########
# HELPERS #
###########

# /*
#  * Returns a line from a file connection and stores in a pocket S4 object
#  * the gene name, chromosome name start and end coordinates as well as the
#  * strand orientation. This function is implemented in Reload.
#  */
getNextLine <- function(con){
	line <- readLines(con, n=1)
	obj <- "NA"
	if (length(line) > 0) {
		line <- strsplit(line, split="\t")[[1]]
		geneName <- line[9]
		geneName <- strsplit(geneName, split=";")[[1]]
		geneName <- geneName[grepl("Name=", geneName)]
		geneName <- strsplit(geneName, "=")[[1]][2]
		obj <- NEW.POCKET(line[1], name=geneName, as.numeric(line[4]), as.numeric(line[5]), line[7])
	}
	obj
}

# /* 
#  * Reloads the list of pockets of coordinates data with a new line.
#  * The data stored in c is passed to b, and data from b is passed to a
#  * and a new line of coordinates data is stored at c.
#  */
Reload <- function(listOfPockets=NULL, con){
	lop <- listOfPockets
	if (is.null(lop)) {
		lop <- new("listOfPockets", list=list())
	}
	if (length(lop@list) == 0) {
		lop@list[["c"]] <- getNextLine(con)
		return(lop)
	}
	if (length(lop@list) == 1) {
		lop@list[["b"]] <- lop@list[["c"]]
		lop@list[["c"]] <- getNextLine(con)
		return(lop)
	}
	if (length(lop@list) >= 2) {
		lop@list[["a"]] <- lop@list[["b"]]
		lop@list[["b"]] <- lop@list[["c"]]
		lop@list[["c"]] <- getNextLine(con)
		return(lop)
	}
}


getGFF <- function(lopos) {
	gff <- c()
	for (name in names(lopos@list)) {
		entry <- lopos@list[[name]]
		comments <- paste("Name=", entry@name, sep="")
		gff <- rbind(gff, c(entry@chrname, ".", "Promoter", entry@coords[1], entry@coords[2], ".", entry@strand, ".", comments))
	}
	gff
}

write.gff <- function(lopos, file) {
	print("Formatting data.")
	gff <- getGFF(lopos)
	print("Writing gff object to file.")
	write.table(gff, file, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
	print("Finished!")
}

write.fasta <- function(lopos, file) {
	print("Collapsing lits of pockets of sequences object.")
	lopos <- collapse.lopos(lopos)
	print("Writing DNAStringSet object to file in fasta format.")
	writeXStringSet(lopos, filepath=file, format="fasta")
	print("Finished!")
}

collapse.lopos <- function(lop) {
	if(class(lop@list[[1]]) != "pocketOfSeq"){
		stop("lop object does not contain pocketOfSeq")
	}
	gnames <- names(lop@list)
	min <- 0
	max <- length(gnames)
	pb <- txtProgressBar(min, max, style=3)
	sequences <- DNAStringSet()
	for (name in names(lop@list)) {
		sequences <- c(sequences, lop@list[[name]]@seq)
		min <- min + 1
		setTxtProgressBar(pb, min)
	}
	names(sequences) <- gnames
	sequences
}

# /*
#  * Could this script be more efficient
#  */
getPromSeqs <- function(lop, sequences, max_size, begin=FALSE) {
	los <- new("listOfPockets", list=list())
	b.strand <- lop@list[["b"]]@strand
	b.start <- lop@list[["b"]]@start
	b.end <- lop@list[["b"]]@end
	b.name <- lop@list[["b"]]@name
	b.chrname <- lop@list[["b"]]@chrname
	c.start <- lop@list[["c"]]@start
	chr.start <- 1
	chr.end <- width(sequences[names(sequences) %in% b.chrname])

	if (b.start <= 0 | b.end <= 0 | c.start <= 0) {
		return(NULL)
	}

	if (b.start > chr.end) {
		return(NULL)
	}
	if (b.end > chr.end) {
		return(NULL)
	}

	if (begin) {
		if (b.strand == "+"){
			if (b.start - chr.start > max_size) {
				pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = subseq(sequences[names(sequences) %in% b.chrname], b.start - max_size, b.start-1), coords = c(b.start - max_size, b.start-1), strand = "+")
				return(pos)
			}
			else {
				pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = subseq(sequences[names(sequences) %in% b.chrname], chr.start, b.start-1), coords = c(chr.start, b.start-1), strand = "+")
				return(pos)
			}
		}
		# IF b.strand != "-"
		else {
			if (c.start - b.end > max_size | c.start - b.end < 0) {
				if (b.end +  max_size <= chr.end) {
					pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = reverseComplement(subseq(sequences[names(sequences) %in% b.chrname], b.end, b.end + max_size)), coords = c(b.end, b.end + max_size), strand = "-")
					return(pos)
				}
				else {
					pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = reverseComplement(subseq(sequences[names(sequences) %in% b.chrname], b.end+1, chr.end)), coords = c(b.end+1, chr.end), strand = "-")
					return(pos)
				}
			}
			else {
				pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = reverseComplement(subseq(sequences[names(sequences) %in% b.chrname], b.end, c.start)), coords = c(b.end, c.start), strand = "-")
				return(pos)
			}
		}
	}

	# IF BEGIN = FALSE
	else {
		a.end <- lop@list[["a"]]@end
		if (b.start <= 0 | b.end <= 0 | c.start <= 0 | a.end <= 0) {
			return(NULL)
		}
		if (b.strand == "+"){
			if (b.start - a.end > max_size | b.start - a.end < 0) {
				if (b.start - max_size > 0) {
					pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = subseq(sequences[names(sequences) %in% b.chrname], b.start - max_size, b.start-1), coords = c(b.start - max_size, b.start-1), strand="+")
					return(pos)
				}
				else {
					pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = subseq(sequences[names(sequences) %in% b.chrname], chr.start, b.start-1), coords = c(chr.start, b.start-1), strand = "+")
					return(pos)
				}
			}
			else {
				pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = subseq(sequences[names(sequences) %in% b.chrname], a.end, b.start-1), coords = c(a.end, b.start-1), strand = "+")
				return(pos)
			}
		}
		# IF b.strand != "-"
		else {
			if (c.start - b.end > max_size | c.start - b.end < 0) {
				if (b.end + max_size <= chr.end) {
					pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = reverseComplement(subseq(sequences[names(sequences) %in% b.chrname], b.end+1, b.end + max_size)), coords = c(b.end+1, b.end + max_size), strand = "-")
					return(pos)
				}
				else {
					pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = reverseComplement(subseq(sequences[names(sequences) %in% b.chrname], b.end+1, chr.end)), coords = c(b.end+1, chr.end), strand = "-")
					return(pos)
				}
			}
			else {
				if (c.start <= chr.end) {
					pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = reverseComplement(subseq(sequences[names(sequences) %in% b.chrname], b.end+1, c.start)), coords = c(b.end+1, c.start), strand = "-")
					return(pos)
				}
				else {
					pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = reverseComplement(subseq(sequences[names(sequences) %in% b.chrname], b.end+1, chr.end)), coords = c(b.end+1, chr.end), strand = "-")
					return(pos)
				}
			}
		}
	}
}

finish_it <- function(lop, sequences, max_size) {
	print("finish it! Cano wins. tuc tuc, tuc tuc, tuc tuc.")
	NULL
}


###########
# METHODS #
###########

setMethod("show", signature("listOfPockets"), function(lop) {
	for (name in names(lop@list)) {
		show(lop@list[[name]])
	}
})

setMethod("show", signature("pocket"),
function(p) {
	str <- paste(names(p), p@chrname, p@name, p@start, p@end, p@strand, sep="   ")
	print(str)
})

setMethod("show", signature("pocketOfSeq"),
function(p) {
	str <- paste(p@name, p@chrname, sep="|")
	print(str)
	print("Coordinates")
	print(p@coords)
	print("Strand")
	print(p@strand)
	print(p@seq)
	print("------------##  CAT  ##--------------")
})


##################
# MAIN FUNCTIONS #
##################

extractPromoters <- function(input, sequences, max_size){
	con <- file(input, "r")
	lenfile <- length(readLines(con))
	close(con)
	con <- file(input, "r")
	if (!"file" %in% class(con) | class(sequences) != "DNAStringSet" | class(max_size) != "numeric") {
		stop("Input data are corrupted.\nInput data must be file connection and DNAStringSet objects\n")
	}
	lop <- Reload(con=con)
	los <- new("listOfPockets", list=list())
	min <- 0
	max <- lenfile
	pb <- txtProgressBar(min, max, style=3)
	while (TRUE) {
		if ("c" %in% names(lop@list)) {
			if (is.atomic(lop@list[["c"]])) {
				if (lop@list[["c"]] == "NA") {
					gname <- lop@list[["b"]]@name
					los@list[[gname]] <- finish_it(lop, sequences, max_size)
					close(con)
					return(los)
				}
			}
			if ("b" %in% names(lop@list)) {
				gname <- lop@list[["b"]]@name
				if ("a" %in% names(lop@list)) {
					los@list[[gname]] <- getPromSeqs(lop, sequences, max_size)
#					show(lop)
					lop <- Reload(lop, con=con)
				}
				else {
					los@list[[gname]] <- getPromSeqs(lop, sequences, max_size, begin=TRUE)
#					show(lop)
					lop <- Reload(lop, con=con)
				}
			}
			else {
				lop <- Reload(lop, con=con)
			}
		}
		else {
			lop <- Reload(lop, con)
		}
		min <- min + 1
		setTxtProgressBar(pb, min)
	}
}









