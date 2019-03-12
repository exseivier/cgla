#!/usr/bin/env bash

usage(){
cat << USAGE
	****************************************
	*            scripts folder            *
	*       cglabrata parent folder        *
	*        Projects parent folder        *
	****************************************

	This folder holds the scripts used to process the rna-seq
	data from experiments in "Candida glabrata" experiments
	where this yeast was challenged to grow in stressing media
	(1) nitrosative media and (2) shift pH condition.

	List of scripts:
	---------------

	(1) matchPattern-2.py		This script uses the sequences of a fasta file
					to search the ocurrences and the location of
					a pattern matching on every sequence.
					Every sequence is divided into 5 categories.
					(1) very near, (2) near, (3) moderate, (4) far,
					and (5) very far, and every sub-sequence has the
					same size. Then for every sub-sequence and for every
					gene sequence (promoter, UTR, CDS) this script returns
					the occurrence and the location of the pattern
					matching for every gene and for every subsequence.

					output file format

					gene_name	very_near	near	moderate	far	very_far
					gene1		1			0		3			0
					gene2		0			5		0			0
					...			...			...		...			...
					geneN		0			3		3			3

	(2) cumfreq-2.R			This script takes the output file of matchPattern-2.py and
					plots the cummulative frequency of the logFC for the entire
					dataset and for the sub-datasets which memebers have at least
					one pattern match in one of all subsequences (very_near, near,
					moderate, far and very_far). Also it calculates the wilcoxon
					rank-sum test to evaluate if the shift of the means of the
					subpopulations (near, far, very_far etc...) compared to the
					entire population is significative.

	(3) gffPromCreat.sh		This script takes a gff file and creates the coordinates
					of the intergenic regions based on "gene" features.
	
	(4) pickRandomSeqs.sh		This script takes a fasta file name and an integer
					number "N" as arguments, and returns a fasta file with "N" randomly
					selected sequences.
	
	(5) obo-tools.py		This script takes a GOIDs list and collapses them into 
					a more general (parental) GO terms list


USAGE
}

usage
