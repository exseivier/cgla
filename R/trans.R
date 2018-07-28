
#############
# LIBRARIES #
#############

# // NO LIBRARIES
# // LIBRARIES LOADED FROM extract.R

###########
# OBJECTS #
###########

# // NO OBJECTS

###########
# HELPERS #
###########

# /* 
#  * Takes the *.tab files of annotations from YGOB database
#  * www.ygob.ucd.ie/ygob/data/v7-Aug2012
#  * and returns a gff formated file
#  */
ygobTab2gff <- function(ygobTab_file, faFile) {
	table <- read.delim(ygobTab_file, header=FALSE, as.is=TRUE, sep="\t", quote="")
	strand <- ifelse(table$V2 == 1, "+", "-")
	table$V2 <- strand
	genome <- readDNAStringSet(faFile, format="fasta")
	
	# TESTING INSECTICIDE for BUG No. 1
	table$V6 <- unlist(lapply(table$V6, function(x) names(genome[grepl(paste("_", x, "$", sep=""), names(genome))])))
	
	# BUG No. 1: 
	#table$V6 <- gsub(" .*$", "", names(genome)[table$V6])
	
	len <- length(table$V6)
	comments <- paste("Name=", table$V1, ";SHORTNAME=", table$V7, sep="")
		result <- data.frame(chr=table$V6, source=rep(".", len), type=rep("mRNA", len), start=table$V3, end=table$V4, unk1=rep(".", len), strand=table$V2, phase=rep(".", len), comments=comments)
		out <- strsplit(ygobTab_file, "\\.")[[1]]
	out <- paste(out[1:(length(out)-1)], collapse=".")
	out <- paste(out, ".gff", sep="")
	write.table(result, file=out, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
	print("Data was stored in gff format")
}

###########
# METHODS #
###########

# // NO METHODS
