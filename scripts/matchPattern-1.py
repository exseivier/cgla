#!/usr/bin/env python

from sys import argv, exit

def DNASeqSet(filename):
	"""(STR) -> HASH[STR:STR]
	From a fasta file it takes every sequence and stores them into a hash data structure
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

def revComp(dna):
	"""(STR) -> STR
	Returns the reverse complementary sequence of a DNA
	"""
	rc = ""
	nuc_dict = {
				"A":"T",
				"C":"G",
				"G":"C",
				"T":"A",
				}
	for nt in dna:
		rc = nuc_dict[nt] + rc
	
	return rc

def matchStr(pattern, dna, mismatch_allowed):
	"""(STR, HASH[STR:STR]) -> ARRAY[INT]
	Returns an array with the start position where a match of "pattern" ocurred in "dna".
	"""
	i = 0
	array = []
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
	filename = argv[1]
	output = ".".join(filename.split(".")[0:-1]) + ".cat"
	pattern = argv[2]
	mismatch_allowed = int(argv[3])
	output = pattern + "_" + str(mismatch_allowed) + "m.cat"
	print "Uploading sequences..."
	DNA = DNASeqSet(filename)
	
	print "Computing matches..."
	matches = {}
	for header, sequence in DNA.iteritems():
		positions = matchStr(pattern, sequence, mismatch_allowed)
		positions_rc = matchStr(revComp(pattern), sequence, mismatch_allowed)
		positions.extend(positions_rc)
		matches[header] = positions
	
	print "Categorising..."
	categories = {}
	FO = open(output, "w+")
	FO.write("gene_name\tvery_near\tnear\tmoderate\tfar\tvery_far\n")
	for header, positions in matches.iteritems():
		categories[header] = {
								"vn":0,
								"n":0,
								"m":0,
								"f":0,
								"vf":0,
								}
		print positions
		for position in positions:
			if position > 0 and position < 300:
				print "We got a less 300"
				categories[header]["vn"] += 1	
			elif position >= 300 and position < 600:
				categories[header]["n"] += 1
			elif position >= 600 and position < 900:
				categories[header]["m"] += 1
			elif position >= 900 and position < 1200:
				categories[header]["f"] += 1
			elif position >= 1200 and position <= 1501 - len(pattern):
				categories[header]["vf"] += 1
			else:
				print "ERROR! Position ", position, " was not categorized!"
				#exit(1)

		FO.write("%s\t%d\t%d\t%d\t%d\t%d\n" % (header, categories[header]["vn"], categories[header]["n"], categories[header]["m"], categories[header]["f"], categories[header]["vf"]))
	return 0

if __name__ == "__main__": main()



