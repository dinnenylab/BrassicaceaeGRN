#!/usr/bin/env python
synopsis = "\n### usage: find_elements_in_promoters.py <cis-elements.list> <promoters.fasta> <output.txt>\n\
### find locations of elements in <cis-elements.list>, within each sequence in <promoters.fasta>, and print to <output.txt>\n\
### <cis-elements.list> contains an element ID and an element sequence, delimited by a tab, per each line \n\
### element sequences can be either sequences or python regular expression patterns \n\
### <promoters.fasta> should be formatted as one-line sequences (e.g. using ""fasta-formatter"") \n\
### <output.txt> contains the promoter ID, start and end positions, the matched sequence, the element ID, and the element pattern, ... \n\
### ... tab-delimited and one occurrence per line \n\
### searches only for the top strand. add inverted sequences to the <cis-elements.list> for all elements \n\
### for overlapping matches of the same element, only the first occurrence will be reported. all non-overlapping matches will be reported \n\
### 201114, script renamed to 'find_motifs_in_promoters.py';\n\
### copyleft by ohdongha@gmail.com 150927 ver 2.1, now can accept non-capital sequences;\n\
### copyleft by ohdongha@gmail.com 150927 ver 2.1, now print start and end positions;\n\
### copyleft by ohdongha@gmail.com 150503 ver 2.0, adding capability to deal with regular expression;\n\
### copyleft by ohdongha@gmail.com 150328 ver 1.0\n"

import sys
import re

### helper function
def find_all_pattern(a_str, sub):
	start = 0
	match = ""
	while True:
		match = re.search(sub, a_str[start:])
		if match :
			start = a_str.find(match.group(), start)
			yield (start + 1, start + len(match.group()), match.group())
			start += len(match.group()) # use start += 1 to find overlapping matches
		else : return

### reading in arguments
try: 
	fin_list = open(sys.argv[1], "rU")
	fin_fasta = open(sys.argv[2], "rU")
	fout = open(sys.argv[3], "w")
except ValueError :
	print synopsis
	sys.exit(0)
except IndexError :
	print synopsis
	sys.exit(0)


### reading in <cis-elements.list>
lines_in_elementsList = 0
faultyLines_in_elementsList = 0
elementID = ''
element_dict = dict()

print "\nReading ", sys.argv[1].strip(), ": \n"
for line in fin_list:
	tok = line.split('\t')
	lines_in_elementsList = lines_in_elementsList + 1
	try:
		if (len(tok[1].strip()) > 3):
			elementID = tok[0].strip()
			element_dict[elementID] = tok[1].strip().upper()
		else:
			faultyLines_in_elementsList = faultyLines_in_elementsList +1			
	except IndexError :
		faultyLines_in_elementsList = faultyLines_in_elementsList +1
	except ValueError :
		faultyLines_in_elementsList = faultyLines_in_elementsList +1

print "Out of total ", lines_in_elementsList, " lines in ", sys.argv[1].strip(), ", ", faultyLines_in_elementsList, " lines were rejected."
fin_list.close()


### reading in <promoters.fasta>, counting, and printing
promoterID = ""
promoterSeq = ""
newSeq = 0
ElementsFound = []
counter = 0

print "\nReading ", sys.argv[2].strip(), ": \n"
for line in fin_fasta:
	if (line[0] == '>'):
		promoterID = line[1:].strip()
		newSeq = 1
	elif ( newSeq == 1):
		promoterSeq = line.strip().upper()
		for key in sorted(element_dict):
			ElementsFound = list(find_all_pattern(promoterSeq, element_dict[key]))
			for item in ElementsFound :
				fout.write(promoterID + '\t' \
						+ str(item[0]) + '\t' \
						+ str(item[1]) + '\t' \
						+ str(item[2]) + '\t' \
						+ key + '\t' \
						+ element_dict[key] + '\n' )
			#print "\nElement %s was processed. \n" % ( key ) 
		newSeq = 0 
		counter = counter + 1
		if ( counter % 1000 == 0):
			sys.stdout.write("\r   processing ~%d promoters" % int(counter))
			sys.stdout.flush()
	
print "\ndone"
fin_fasta.close()
fout.close()
