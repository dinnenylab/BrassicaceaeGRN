#!/usr/bin/env python
import sys, re

synopsis = "\n\
find_motifs_in_promoters.py <cis-elements.list> <promoters.fa> <output.txt>\n\
 - find locations of elements in <cis-elements.list>, within each sequence in\n\
   <promoters.fa>, and print to tab-delimited <output.txt>;\n\
 - <cis-elements.list> contains an element ID and an element pattern, tab-\n\
   delimited, one element per line; an element pattern can be either a sequence\n\
   or a python regular expression pattern;\n\
 - <promoters.fa> should be formatted as one-liner sequences (e.g. using\n\
   'fasta-formatter'); in a sequence header (ID), only the first field before\n\
   the first space or tab is used as the promoter ID;\n\
 - <output.txt> contains the promoter ID, followed by start and end positions,\n\
   the matched sequence, the element ID, and the element pattern, tab-delimited\n\
   and one match per line \n\
 - searches only the top strand of <promoter.fa>, hence add inverted sequences\n\
   to <cis-elements.list> for all elements;\n\
 - for overlapping matches of the same element, only the first occurrence are\n\
   reported, while all non-overlapping matches are reported;\n\
by ohdongha@gmail.com ver2.2.1 20210121\n"
#version_history
#210121 ver 2.2.1 # minimally modified to work with python 3
#201226 ver 2.2 # use only the first field of sequence IDs in <promoter.fasta>
#150927 ver 2.1 # now can accept non-capital sequences
#150927 ver 2.1 # now print start and end positions
#150503 ver 2.0 # adding capability to deal with regular expression
#150328 ver 1.0

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
	fin_list = open(sys.argv[1], "r")
	fin_fasta = open(sys.argv[2], "r")
	fout = open(sys.argv[3], "w")
except:
	print( synopsis )
	sys.exit(0)

### reading in <cis-elements.list>
lines_in_elementsList = 0
faultyLines_in_elementsList = 0
elementID = ''
element_dict = dict()

print( "\nReading ", sys.argv[1].strip(), ": \n" )
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

print( "Out of total ", lines_in_elementsList, " lines in ", sys.argv[1].strip(), ", ", faultyLines_in_elementsList, " lines were rejected." )
fin_list.close()


### reading in <promoters.fasta>, counting, and printing
promoterID = ""
promoterSeq = ""
newSeq = 0
ElementsFound = []
counter = 0

print( "\nReading ", sys.argv[2].strip(), ": \n" )
for line in fin_fasta:
	if (line[0] == '>'):
#		promoterID = line[1:].strip()
		promoterID = line[1:].strip().split()[0] # use only the first field
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
	
print( "\ndone" )
fin_fasta.close()
fout.close()