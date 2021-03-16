#!/usr/bin/env python
import sys, re, argparse
from argparse import RawTextHelpFormatter
	
###################################################
### 0. script description and parsing arguments ###
###################################################
synopsis1 = "\
  accept coordinates of two genomic regions (chrID, start, and end) and mark\n\
  regions in <target> that are included in any region in <query>."

synopsis2 = "detailed description:\n\
 1. Input:\n\
  - <query>: tab-delimited with an unique ID (uID), chrID, start, and end for\n\
     each genomic region; genomic regions can overlap, but uID must be unique;\n\
     header not required\n\
  - <target>: tab-delimited with an chrID, start, and end for each region in\n\
     columns <N>, <N+1>, and <N+2>, respectively; e.g. if the file starts with\n\
     chrID, <N> == 1\n\
  - recomended uID format for <query> is 'chrID:start-end__annotation'; uID is\n\
     the only information appear in the output\n\
  - in general, <query> should have the larger genomic regions than <target>;\n\
     if the other way around, use '-r' option below\n\
 2. Output and options:\n\
  - <output> contains the entire line of <target> and uID, tab-delimited; if\n\
     there is no region in <query> that includes the genomic region, '_na_'\n\
     will be printed instead of uID\n\
  - '-r'|'--reverse': look for regions in <query> INCLUDED within a region in\n\
     <target>; if the target region includes multiple <query> regions, print\n\
     the line in <target> multiple times with each uID; [False]\n\
  - '-c'|'--counter': report progress per every X lines; [10000]\n\
by ohdongha@gmail.com 20200127 ver 1.2\n\n"

#version_history
#20210314 script renamed to "genomic_regions_mark_overlaps.py"
#20201114 script renamed to "genomic_regions_mark_regions_included_in_others.py"
#20200127 ver 1.2 # '-r' option added
#20200125 ver 1.1 # minor modification to the output format
#20150419 ver 1.0 

parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)

# positional parameters
parser.add_argument('query', type=str, help="See below")
parser.add_argument('N', type=int)
parser.add_argument('target', type=str, help="See below")
parser.add_argument('output', type=str, help="See below")

# options
parser.add_argument('-r', '--reverse', action="store_true", default=False)
parser.add_argument('-c', '--counter', type=int, default=10000)

args = parser.parse_args()

fin_query = open(args.query, "rU")
chrID_colIndex = args.N
fin_target = open(args.target, "rU")
fout = open(args.output, "w")


##########################################################
### 1. read <query> and <target> and write to <output> ###
##########################################################
### reading in <query>
lines_in_region_table = 0
faultyLines_in_region_table = 0
uID = ""
chrID = ""
start_q = 0
end_q = 0
regionNum = 0

uID_dict = dict()
chrID_dict = dict()
start_dict = dict()
end_dict = dict()

for line in fin_query:
	tok = re.split('\t', line)
#	tok = line.split('\t')
	lines_in_region_table = lines_in_region_table + 1
	try:
		uID = tok[0].strip()
		chrID = tok[1].strip()	
		start_q = int(tok[2].strip())
		end_q = int(tok[3].strip())
		if end_q >= start_q :
			uID_dict[regionNum] = uID
			chrID_dict[regionNum] = chrID
			start_dict[regionNum] = start_q
			end_dict[regionNum] = end_q
			regionNum = regionNum + 1
		else:
			faultyLines_in_region_table = faultyLines_in_region_table +1			
	except (IndexError, ValueError) :
		faultyLines_in_region_table = faultyLines_in_region_table +1

print "Reading <query> file: %s \n" % fin_query.name
print "Out of total %d line in query, %d were rejected.\n" % (lines_in_region_table, faultyLines_in_region_table)
print "Now marking regions in <target> file: %s \n" % fin_target.name
fin_query.close()


### reding lines in <target> and print to <output>
Lines = 0
collectedLines = 0
chrID = ""
start_t = 0
end_t = 0
counted = 0

for line in fin_target:
	try :
		tok = re.split('\t', line)
		chrID = tok[chrID_colIndex - 1].strip()
		start_t = int( tok[chrID_colIndex].strip() )
		end_t = int( tok[chrID_colIndex + 1].strip() )
		if args.reverse:
			for i in range(0, regionNum) :
				if  (chrID_dict[i] == chrID) and (start_dict[i] >= start_t) and (end_dict[i] <= end_t) : 
					fout.write(line.strip() + '\t' + uID_dict[i] + '\n')
					counted = 1		
		else:
			for i in range(0, regionNum) :
				if  (chrID_dict[i] == chrID) and (start_dict[i] <= start_t) and (end_dict[i] >= end_t) : 
					fout.write(line.strip() + '\t' + uID_dict[i] + '\n')
					counted = 1
		if (counted == 1):
			collectedLines = collectedLines + 1
			counted = 0
		else:
			fout.write(line.strip() + '\t_na_\n') # if not included in (or including) any region in <query>
		Lines = Lines + 1
		if ( Lines % args.counter == 0):
			sys.stdout.write("\r   marked %d regions out of %d" % (collectedLines, Lines))
			sys.stdout.flush()
	except (ValueError, IndexError):
		pass

print "\n\nFor total ", Lines, " regions in ", sys.argv[3].strip(), ", ", collectedLines, " were marked."
print "Printing to ", sys.argv[4].strip(), ": \n"
fin_target.close()
fout.close()

print "done"
