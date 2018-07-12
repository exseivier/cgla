# /* 
#    This methods container are a collection of functions intended
#    to extract feature elements defined in an annotation file
# */ 


# LIBRARIES
library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(BSgenome)

# OBJECTS

# HELPERS

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

# METHODS

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


