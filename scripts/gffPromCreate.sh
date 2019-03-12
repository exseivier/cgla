#!/usr/bin/env bash

#	VARIABLES
GFF=$1	# gff file name
SIZE=$2

#	AWK SCRIPT
awk "
	BEGIN 	{
			bs = \"\";
			ss = \"\";
			as = \"\";
			bstart = 0;
			bend = 0;
			start = 0;
			end = 0;
			astart = 0;
			aend = 0;
			counter = 0;
			}
	{
		if(\$3 == \"mRNA\") {
			bstart = start;
			bend = end;
			bs = ss;
			start = astart;
			end = aend;
			ss = as;
			astart = \$4;
			aend = \$5;
			as = \$0
			counter = counter + 1
		#	print counter
			if(counter >= 2){
				split(ss, ssa, \"\t\");
				split(bs, bsa, \"\t\");
				split(as, asa, \"\t\");
				if (bsa[1] != ssa[1]){
					bstart = 0;
					bend = 0;
				}
		#		print \"---------------\"
		#		print bstart, \"-\", bend
		#		print start, \"-\", end
		#		print astart, \"-\", aend
				split(ssa[9], comments, \";\")
				split(comments[2], head, \"=\")
				
				
				if (ssa[7] == \"+\") {
					if (bsa[1] != ssa[1]) {
						bstart = 0;
						bend = 0;
					}
					if (bend + 1 > start - 1) {
						if (start - $2 > 0) {
							pstart = start - $2
							pend = start - 1
						}
						else {
							pstart = 1
							pend = start - 1
						}
					}
					else {
						pstart = bend + 1
						pend = start - 1
					}

				#	print \"pend - pstart\"
				#	print pend, \"-\" pstart
					if (pend - pstart > $2){
						if(pstart < 0){
							printf \"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n\", ssa[1], head[2], head[2], 1, pend, \".\", \"+\", \".\", ssa[9]
						}
						else {
							printf \"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n\", ssa[1], head[2], head[2], pend - $2, pend, \".\", \"+\", \".\", ssa[9]
						}
					}
					else if (pend - pstart < $3) {
						if(pend - $2 < 0){
							#printf \"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n\", ssa[1], head[2], head[2], 1, pend, \".\", \"+\", \".\", ssa[9]
						}
						else {
							printf \"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n\", ssa[1], head[2], head[2], pend - $2, pend, \".\", \"+\", \".\", ssa[9]
						}
					}
					else{
						printf \"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n\", ssa[1], head[2], head[2], pstart, pend, \".\", \"+\", \".\", ssa[9]
					}
#					printf \"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n\", ssa[1], head[2], ssa[3], ssa[4], ssa[5], \".\", ssa[7], \".\", ssa[9]
				}


				else if (ssa[7] == \"-\") {
					pstart = end + 1
					if (asa[1] != ssa[1]){
						pend = pstart + $2 + 1
					}
					if (end +1 > astart -1) {
						pend = end + $2
					}
					else {
						pend = astart - 1
					}
				#	print \"pend - pstart\"
				#	print pend, \"-\", pstart
#					printf \"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n\", ssa[1], head[2], ssa[3], ssa[4], ssa[5], \".\", ssa[7], \".\", ssa[9]
					if (pend - pstart > $2){
						printf \"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n\", ssa[1], head[2], head[2], pstart, pstart + $2, \".\", \"-\", \".\", ssa[9]
					}
					else if (pend - pstart < $3) {
						printf \"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n\", ssa[1], head[2], head[2], pstart, pstart + $2, \".\", \"-\", \".\", ssa[9]
					}
					else {
						printf \"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n\", ssa[1], head[2], head[2], pstart, pend, \"this\", \"-\", \".\", ssa[9]
					}
				}
			}
		}
	}
" $GFF
