# /* 
#    This methods container is a collection of functions intended
#    to extract feature elements defined in an annotation file
# */ 

#############
# LIBRARIES #
#############

library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(BSgenome)

###########
# OBJECTS #
###########

# // NO OBJECTS

###########
# HELPERS #
###########

# /* 
#  * Using a gff file, it calculates the intergenic coordinates
#  * and returns a gff file with those corrdinates.
#  * It requires a gff file and a fasta file (character vectors)
#  * also it requires a integer vector with the maximum size of
#  * the intergenic region to return.
#  * This function returns a GenomicRanges object
#  */
get_intergenic_coordinates <- function(gff_file, FaFile, max_size) {
	gff <- import.gff(gff_file)
	gff.mRNA <- gff[gff$type == "mRNA"]
	info <- seqinfo(scanFaIndex(FaFile))
	gff.flank <- flank(gff.mRNA, max_size)
	seqlevels(gff.flank) <- seqlevels(info)
	seqinfo(gff.flank) <- info
	gff.flank <- trim(gff.flank, use.names=TRUE)
	gff.dif <- setdiff(gff.flank, gff.mRNA)
	gff.dif.plus <- gff.dif[strand(gff.dif) == "+"]
	gff.dif.minus <- gff.dif[strand(gff.dif) == "-"]
	minus.overlaps <- findOverlaps(flank(gff.dif.minus, 10, start=FALSE), gff.mRNA, select="first")
	plus.overlaps <- findOverlaps(flank(gff.dif.plus, 10, start=FALSE), gff.mRNA, select="first")
	gff.dif.plus$Name <- gff.mRNA[plus.overlaps]$Name
	gff.dif.minus$Name <- gff.mRNA[minus.overlaps]$Name
	gff.dif <- c(gff.dif.plus, gff.dif.minus)
	gff.dif
}

# /* 
#  * Returns a DNAStringSet object with the intergenic regions specified
#  * in a GenomicRanges object (granges). Sequences are extracted from
#  * a DNAStringSet (BSgenome). You could select only a set of genes
#  * to download only this set of genes setting (geneSet).
#  */
get_intergenic_sequences <- function(granges, BSgenome, geneSet){
	if (! class(granges) == "GRanges"){
		stop("ERROR!: GRanges object is corrupted!")
	}
	if (!is.null(geneSet)) {
		granges <- granges[granges$Name %in% geneSet]
	}
	interSeqs <- DNAStringSet()
	pb <- txtProgressBar(0, length(granges), style=3)
	for (i in 1:length(granges)) {
		interSeqs <- c(interSeqs, getSeq(BSgenome, granges[i]))
		setTxtProgressBar(pb, i)
	}
	names(interSeqs) <- granges$Name
	interSeqs
}

# /* 
#  * Returns a relational list object with the names of the
#  * DNAStringSetList object (dssl) and the column positions where
#  * the orthologues names are in orthologues database table
#  * (orthoDB)
#  */
get_positions_table <- function(dssl, orthoDB) {
	positions <- list()
	for(name in names(dssl)){
		for(i in 1:length(orthoDB[1,])){
			if (sum(orthoDB[,i] %in% names(dssl[[name]])) > 0) {
				if (name %in% names(positions)) {
					positions[[name]] <- c(positions[[name]], i)
				}
				else {
					positions[[name]] <- i
				}
			}
		}
	}
	positions
}

# /* 
#  * Using the orthologues database table (orthoDB), this function returns
#  * a DNAStringSetList only with the sequences of the orthologues genes
#  * defined in the orthoogues table (orthoDB).
#  */
get_orthogenes_from_DNAStringSetList <- function(dssl, orthoDB) {
	if(class(dssl) != "list" | class(orthoDB) != "data.frame") {
		stop("Input data tables are corrupted")
	}
	positions <- get_positions_table(dssl, orthoDB)
	print(positions)
	for (name in names(dssl)) {
		gNames <- as.vector(as.matrix(orthoDB[,positions[[name]]]))
		dssl[[name]] <- dssl[[name]][names(dssl[[name]]) %in% gNames]
	}
	dssl
}

# /* 
#  * Takes the directories of the path where the AME output is and load every
#  * ame.tsv file for every analyes performed. AME is a motif enrichment analysis
#  * software which belongs to the MEME suite package.
#  * This function takes a character vector (path) with the path where AME result
#  * directories are and for every directory takes the ame.tsv file and load the
#  * table in a data.frame object.
#  */
load.ame_out <- function(path) {
	ame_out <- list()
	dirs <- list.dirs(path, full.names=FALSE)[list.dirs(path, full.names=FALSE) != ""]
	for (dir in dirs) {
		con <- paste(path, dir, "ame.tsv", sep="/")
		table <- read.delim(con, header=TRUE, as.is=TRUE, sep="\t", quote="", comment.char="#")
		table <- table[order(table[,2]),]
		table$log.p <- -log10(table[,6])
		ame_out[[dir]] <- table
	}
	ame_out
}

###########
# METHODS #
###########

# /* 
#  * This method wraps the get_intergenic_coordinates, get_intergenic_sequences
#  * in order to obtain the intergenic sequences in a DNAStringSet object
#  * or the intergenic regions coordinates in a GenomicRanges object.
#  * This function rquires a gff file, a fasta file, a maximum intergenic size,
#  * and a character vector with the names of the genes you want to select (geneSet).
#  * retSeq is a logic argument: if it is set to TRUE, this function returns the sequences.
#  * If it is turned FALSE, it returns only the intergenic coordinates in a GenomicRanges object.
#  */
setGeneric("get_intergenic_regions", function(gff_file, FaFile, max_size, geneSet=NULL, retSeq=TRUE) standarGeneric("get_intergenic_regions"))
setMethod("get_intergenic_regions", signature("character", "character", "integer"),
	function(gff_file, FaFile, max_size, geneSet=NULL, retSeq=TRUE) {
		cat("\nCalculating intergenic coordinates\n")
		gff.dif <- get_intergenic_coordinates(gff_file, FaFile, max_size)
		if (retSeq) {
			cat("Loading genome\n")
			BSgenome <- readDNAStringSet(FaFile, format="fasta")
			cat("Extracting intergenic sequences\n")
			interSeqs <- get_intergenic_sequences(gff.dif, BSgenome, geneSet)
			cat("\nEverything is OK\n")
			return(interSeqs)
		}
		else {
			return(gff.dif)
		}
})


