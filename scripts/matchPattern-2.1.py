#!/usr/bin/env python

from sys import argv, exit

def DNASeqSet(filename):
	"""(STR) -> HASH[STR:STR]
	From a fasta file, it takes every sequence and stores them into a hash data structure
	"""
	FH = open(filename, "r")
	HASH = {}

	for LINE in FH:
		LINE = LINE.strip()
		if LINE[0] == ">":
			name = LINE
		else:
			try:
				HASH[name] += LINE
			except Exception:
				HASH[name] = LINE
	
	return HASH

def revComp(patterns):
	"""(STR) -> STR
	Returns the reverse complementary sequence of a DNA
	"""
	reversed_sequences = []
	for pattern in patterns:
		rc = ""
		nuc_dict = {
				"A":"T",
				"C":"G",
				"G":"C",
				"T":"A",
				}
		for nt in pattern:
			rc = nuc_dict[nt] + rc

		reversed_sequences.append(rc)
	
	return reversed_sequences

def matchStr(patterns, dna, mismatch_allowed):
	"""(STR, HASH[STR:STR]) -> ARRAY[INT]
	Returns an array with the start position where a match of "pattern" ocurred in "dna".
	"""
	i = 0
	array = []
	for pattern in patterns:
		while i <= len(dna)-len(pattern):
			bad_score = 0
			dna_substr = dna[i:i+len(pattern)]
			for j in xrange(0, len(dna_substr)):
				if dna_substr[j].upper() != pattern[j].upper():
					bad_score += 1
				if bad_score > mismatch_allowed:
					i += 1
					break

			if bad_score <=  mismatch_allowed:
				array.append(int(i+1))
				i += 1
			else:
				i += 1
	
	return array
			


def main():
	"""
	MAIN FUNCTION
	"""
	filename = argv[1] # Fasta file of sequences 
	#output = ".".join(filename.split(".")[0:-1]) + ".cat" # Commented
	patterns = argv[2] # kmer to search for in sequences
	patterns = patterns.split(" ")
	
	#Selecting the max length pattern
	plen = 0 #Pattern length
	mlen = 0 #Max pattern length
	for item in patterns:
		plen = len(item)
		if plen > mlen:
			mlen = plen
		else:
			pass
	
	mismatch_allowed = int(argv[3]) # mismatch allowed
	output = str(patterns[0]) + str(mismatch_allowed) + "m.cat" # output file name
	print "Uploading sequences..."
	DNA = DNASeqSet(filename) # Preparing the HASH with the sequences and headers
	
	print "Computing matches..."
	matches = {}
	for header, sequence in DNA.iteritems():
		positions = matchStr(patterns, sequence, mismatch_allowed)
		positions_rc = matchStr(revComp(patterns), sequence, mismatch_allowed)
		positions.extend(positions_rc)
		matches[header] = {"positions":positions, "lenseq":len(sequence)}
	
	print "Categorising..."
	categories = {}
	FO = open(output, "w+")
	FO.write("gene_name\tvery_near\tnear\tmoderate\tfar\tvery_far\n")
	for header, items in matches.iteritems():
		categories[header] = {
								"vn":0,
								"n":0,
								"m":0,
								"f":0,
								"vf":0,
								}
		positions = items["positions"]
		lenseq = items["lenseq"]
		size = lenseq//5
		p1 = 0 + size
		p2 = p1 + size
		p3 = p2 + size
		p4 = p3 + size
		print positions
		for position in positions:
			if position > 0 and position < p1:
				categories[header]["vn"] += 1	
			elif position >= p1 and position < p2:
				categories[header]["n"] += 1
			elif position >= p2 and position < p3:
				categories[header]["m"] += 1
			elif position >= p3 and position < p4:
				categories[header]["f"] += 1
			elif position >= p4 and position <= lenseq - mlen: # Max length pattern
				categories[header]["vf"] += 1
			else:
				print "ERROR! Position ", position, " was not categorized!"
				#exit(1)

		FO.write("%s\t%d\t%d\t%d\t%d\t%d\n" % (header, categories[header]["vn"], categories[header]["n"], categories[header]["m"], categories[header]["f"], categories[header]["vf"])) # Writig out to output file
	return 0

if __name__ == "__main__": main()



