#!/usr/bin/env python

from sys import exit, argv

debug = False

def create_obo(filename):
	"""(STR) -> HASH
	Creates an OBO object from obo file of Gene Ontology
	"""
	HASH = {}
	goid = ""
	try:
		FH = open(filename, "r")
	except Exception:
		usage()
		print "File ", filename, " does not exist!"
		exit(1)
	
	for line in FH:
		line = line.strip()
		if line.split(" ")[0] == "id:":
			goid = line.split(" ")[1]
			HASH[goid] = {}
		else:
			if goid != "" and line != "[Term]":
				items = line.split(" ")
				if len(items) > 1:
					key = items[0]
					value = " ".join(items[1:])
				else:
					continue
				try:
					HASH[goid][key].append(value)
				except Exception:
					HASH[goid][key] = [value]
			else:
				pass
	FH.close()
	return HASH

def parsing_GOTerm(GOTerm_filename):
	"""(STR) -> ARRAY(CHR)
	Takes the GO Terms from a file and allocates them to an array structure
	"""
	try:
		FH = open(GOTerm_filename, "r")
	except Exception:
		usage()
		print "File ", GOTerm_filename, " does not exist!"
		exit(1)
	GOTerms = []
	for line in FH:
		line = line.strip()
		GOTerms.append(line)
	FH.close()
	return GOTerms


def reduce_obo(OBO_object, GOTerms, namespace):
	"""(HASH[OBO_Object], ARRAY) -> ARRAY
	Returns the parents of the GOTerms
	"""
	reduced_GOID = []
	Red_OBO_object = {}
	for GOID_IN_OBO, items in OBO_object.iteritems():
		if items["namespace:"][0] == namespace:
			Red_OBO_object[GOID_IN_OBO] = items
		else:
			pass
	#END-FOR-LOOP
	
	# We got the reduced OBO object only with the GOIDs from the specified namespace
	# [molecular function, cellular component, biological process]
	red_GOID = []
	for GOID in GOTerms:
		GOID = GOID.strip()
		for GOID_IN_OBO, items in Red_OBO_object.iteritems():
			if GOID == GOID_IN_OBO:
				if "is_a:" in items:
					red_GOID.append(items["is_a:"][:])
				if "part_of:" in items:
					red_GOID.append(items["part_of:"][:])
			else:
				pass
	return red_GOID


def usage():
	"""
		Print usage information.
	"""
	print """
	
	USAGE:
	------
	program <go_categories.obo> <GoTerms_file.txt> <"namespace">
	
	go_categories.obo		All GO categories in OBO format downloaded from Gene Ontology consortium
	GoTerms_file.txt		List of GO terms to filter
	"namespace"			NameSpace to filter All GO categories ["molecular_function" | "cellular process" | "molecular complex"]
	
	"""


def main():
	"""
	MAIN FUNCTION
	"""
	if len(argv) != 4:
		usage()
		exit(1)
	
	OBO_FILE = argv[1]
	GOT_FILE = argv[2]
	NAMESPAC = argv[3]
	OBO = create_obo(OBO_FILE)
	if debug:
		for key, value in OBO.iteritems():
			print key, value
			if "is_a:" in value:
				print "is_a"
				print True
			else:
				print "There is not \"is_a\""
				print False
			if "part_of:" in value:
				print "part_of"
				print True
			else:
				print "There is not a \"part_of\""
				print False
			break
	GOTerms = parsing_GOTerm(GOT_FILE)
	GOTerms_reduced = reduce_obo(OBO, GOTerms, NAMESPAC)
	GOTerms_reduced_uniq = []
	for item in GOTerms_reduced:
		for goterms in item:
			if goterms in GOTerms_reduced_uniq:
				continue
			else:
				GOTerms_reduced_uniq.append(goterms)
	for items in GOTerms_reduced_uniq:
		print items

if __name__ == "__main__": main()
