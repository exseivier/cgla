#############
# LIBRARIES #
#############

library(S4Vectors)
library(Biostrings)

###########
# OBJECTS #
###########

setClass("pocket", representation(chrname="character", name="character", start="numeric", end="numeric", strand="character"))
setClass("pocketOfSeq", representation(chrname="character", name="character", seq="DNAStringSet"))
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

	if (b.start > chr.end) {
		return(NULL)
	}
	if (b.end > chr.end) {
		return(NULL)
	}

	if (begin) {
		if (b.strand == "+"){
			if (b.start - chr.start > max_size) {
				pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = subseq(sequences[names(sequences) %in% b.chrname], b.start - max_size, b.start-1))
				return(pos)
			}
			else {
				pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = subseq(sequences[names(sequences) %in% b.chrname], chr.start, b.start))
				return(pos)
			}
		}
		# IF b.strand != "-"
		else {
			if (c.start - b.end > max_size | c.start - b.end < 0) {
				pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = reverseComplement(subseq(sequences[names(sequences) %in% b.chrname], b.end, b.end + max_size)))
				return(pos)
			}
			else {
				pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = reverseComplement(subseq(sequences[names(sequences) %in% b.chrname], b.end, c.start)))
				return(pos)
			}
		}
	}

	# IF BEGIN = FALSE
	else {
		a.end <- lop@list[["a"]]@end
		if (b.strand == "+"){
#			print(a.end)
#			print(b.start)
#			print(b.start-max_size)
			if (b.start - a.end > max_size | b.start - a.end < 0) {
				if (b.start - max_size > 0) {
					pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = subseq(sequences[names(sequences) %in% b.chrname], b.start - max_size, b.start))
					return(pos)
				}
				else {
					pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = subseq(sequences[names(sequences) %in% b.chrname], chr.start, b.start))
					return(pos)
				}
			}
			else {
				pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = subseq(sequences[names(sequences) %in% b.chrname], a.end, b.start))
				return(pos)
			}
		}
		# IF b.strand != "-"
		else {
			if (c.start - b.end > max_size | c.start - b.end < 0) {
#				print(b.end)
#				print(c.start)
#				print(b.end + max_size)
#				print(chr.end)
				if (b.end + max_size <= chr.end) {
				#	print ("Is here")
					pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = reverseComplement(subseq(sequences[names(sequences) %in% b.chrname], b.end, b.end + max_size)))
					return(pos)
				}
				else {
					pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = reverseComplement(subseq(sequences[names(sequences) %in% b.chrname], b.end, chr.end)))
					return(pos)
				}
			}
			else {
				if (c.start <= chr.end) {
					pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = reverseComplement(subseq(sequences[names(sequences) %in% b.chrname], b.end, c.start)))
					return(pos)
				}
				else {
					pos <- new("pocketOfSeq", chrname = b.chrname, name = b.name, seq = reverseComplement(subseq(sequences[names(sequences) %in% b.chrname], b.end, chr.end)))
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
	print(p@seq)
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









