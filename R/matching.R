
#############
# LIBRARIES #
#############

# // NO LIBRARIES
# // LIBRARIES ARE LOADED FROM extract.R

###########
# OBJECTS #
###########

# // NO OBJECTS

###########
# HELPERS #
###########

# /* 
#  * Finds kmer matches on sequences. It returns a table (matrix) with
#  * the number of matches per kilobase, the p-value based on poisson
#  * distribution, the negative log10 of p-value, and the FDR for each
#  * p-value (p.adjust using BH method).
#  * This function requires two DNAStringSet objects (kmers and Sequences)
#  * and a integer (mismatch).
#  */
findMatches <- function(kmers, sequences, mismatch) {
	if (class(kmers) != "DNAStringSet" | class(sequences) != "DNAStringSet") {
		stop("ERROR!: kmers or sequences objects are/is corrupted!")
	}
	results <- list()
	kmers_names <- names(kmers)
	for (kname in kmers_names) {
		cat("\nCalculating matches for ")
		cat(kname)
		cat("\n")
		seqs_names <- names(sequences)
		observed_matches <- c()
		min <- 0
		max <- length(seqs_names)
		pb <- txtProgressBar(min, max, style=3)
		for(sname in seqs_names) {
			observed_matches <- c(observed_matches, countPattern(kmers[[kname]], sequences[[sname]], max.mismatch=mismatch, fixed=FALSE) + countPattern(kmers[[kname]], reverseComplement(sequences[[sname]]), max.mismatch=mismatch, fixed=FALSE))
			min <- min + 1
			setTxtProgressBar(pb, min)
		}
		names(observed_matches) <- seqs_names
		observed_matches <- observed_matches / (width(sequences)/1000)
		lam <- mean(observed_matches)
		print(lam)
		observed_matches <- cbind(observed_matches, ppois(lambda=lam, q=observed_matches, log.p=FALSE, lower.tail=FALSE))
		observed_matches <- cbind(observed_matches, -ppois(lambda=lam, q=observed_matches[,1], log.p=TRUE, lower.tail=FALSE)/log(10))
		observed_matches <- cbind(observed_matches, p.adjust(observed_matches[,2], method="BH", n=length(observed_matches[,2])))
		colnames(observed_matches) <- c("matches/Kb", "PVal", "log10_PVal", "FDR")
		results[[kname]] <- observed_matches
	}
	results
}

# /*
#  * Harvets the negative log10 of pvalues for each kmer and for each set of sequences
#  * and returns a matrix with the gene names in rownames and the names of the yeast
#  * species as column names.
#  * It requires a list object with the matches (produced by findMatches function), and
#  * orthologous relational table (data.frame) with the name of the genes an their orthologues
#  * and a positions_table (matrix) with the name of the names(matches) in rownames pointing out
#  * the column positions of the orthologues gene names for this particular species.
#  */
matches2heatmap <- function(matches, orthoDB, position_table) {
	if(class(matches) != "list" | class(orthoDB) != "data.frame" | class(position_table) != "matrix") {
		stop("Input data tables are corrupted")
	}
	positions <- colnames(orthoDB)
	if (is.null(positions)) {
		stop("orthoDB data frame has no column names")
	}
	pvalues_table <- c()
	row_names <- c()
	for (row in 1:length(orthoDB[,1])) {
		pvalues <- c()
		row_names <- c(row_names, orthoDB[row,1])
		for(position in positions) {
			matches_name <- position_table[position_table[,2] %in% position | position_table[,3] %in% position, 1]
			kmer_names <- names(matches[[matches_name]])
			summing <- 0
			for (kmer_name in kmer_names) {
				if (length(matches[[matches_name]][[kmer_name]][rownames(matches[[matches_name]][[kmer_name]]) %in% orthoDB[row, position],2]) == 0) {
					summing <- NA
				}
				else {
					summing <- summing + matches[[matches_name]][[kmer_name]][rownames(matches[[matches_name]][[kmer_name]]) %in% orthoDB[row, position],3]
				}
			}
			pvalues <- c(pvalues, summing)
		}
		pvalues_table <- rbind(pvalues_table, pvalues)
	}
	rownames(pvalues_table) <- row_names
	pvalues_table
}

###########
# METHODS #
###########

# /*
#  *
#  */
setGeneric("match.kmers", function(kmers, sequences) standarGeneric("match.kmers"))
setMethod("match.kmers", signature("DNAStringSet", "DNAStringSet"),
	function(kmers, sequences) {
	matches.table <- findMatches(kmers, sequences)
	# // TODO: To make a wraper to find kmer matches for those functions
})


