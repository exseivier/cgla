# LIBRARIES


# OBJECTS

# HELPERS

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

# METHODS

setGeneric("match.kmers", function(kmers, sequences) standarGeneric("match.kmers"))
setMethod("match.kmers", signature("DNAStringSet", "DNAStringSet"),
	function(kmers, sequences) {
	matches.table <- findMatches(kmers, sequences)

})
