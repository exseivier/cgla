#!/usr/bin/env python

from sys import argv, exit

def main():
	"""

	"""
	filename = argv[1]
	outfile = argv[2]
	FI = open(filename, "r")
	FO = open(outfile, "w+")
	str_out = ""
	for line in FI:
		line = line.strip()
		gene = line.split("\t")[2]
		gots = line.split("\t")[3]
		# Check point 1... [OK]
		#print gene
		#print gots
		gots = gots.split(",")
		for got in gots:
			str_out += "%s\t%s\n" % (gene, got)
	FO.write(str_out)

	FI.close()
	FO.close()

if __name__ == "__main__": main()

