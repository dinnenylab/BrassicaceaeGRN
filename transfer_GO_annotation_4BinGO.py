#!/usr/bin/env python
import sys, argparse
from argparse import RawTextHelpFormatter

###################################################
### 0. script description and parsing arguments ###
###################################################
synopsis1 = "\
 create a custom GO annotation for BinGO, by transfering annotation\n\
 from a reference genome (e.g. Arabidopsis) to an annotated list of genes.\n\
"
synopsis2 = "detailed description:\n\
1. Input, output files, and parameters\n\
 - <anchor_list> : geneID and corresponding refID (anchorID), tab-delimited,\n\
 - <GO_anno> : annotation from GO consortium, e.g. 'gene_association.tair' for\n\
    Arabidopsis thaliana; http://www.geneontology.org/page/download-annotations\n\
 - <output> : custum annotation, as specified by BinGO;\n\
 - '-i Col4anchorIDs': set the column index to locate anchorIDs in the <GO_anno>\n\
    file; default=10; if the <GO_anno> file is in recent .gaf format, you may want\n\
    to use '-i 2' option\n\
2. GO annotation transfer process\n\
 - Assume A. thaliana (At) is the reference GO annotation,\n\
 - <anchor_list> includes the pairing of each geneID to an At locus ID\n\
    (e.g. AT1G10000), associated with the geneID,\n\
 - <anchor_list> can be geneIDs and their At BLAST best hits, for example,\n\
 - from the <GO_anno> file (i.e. 'gene_association.tair'), At locus IDs and\n\
    all GO terms associated with them are extracted,\n\
 - to each geneID, all GO terms annotated for its 'anchor' refID \n\
    (i.e. At locus ID) are transfered in the BinGO annotation format\n\
    (see: https://www.psb.ugent.be/cbd/papers/BiNGO/Customize.html).\n\n\
ohdongha@gmail.com 20210315 ver 0.2\n"

#version_history
#20210315 ver 0.2 mode compatible with python 3
#20180618 ver 0.1 added '-i' option to set the column number to indicate the location of anchorIDs in <GO_anno>
#20170717 ver 0.0

# arguments and parameters
parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)

parser.add_argument('anchor_list', type=argparse.FileType('r'), help="geneID-refID association file, see below")
parser.add_argument('GO_anno', type=argparse.FileType('r'), help="GO annotation file for refIDs, see below")
parser.add_argument('output', type=argparse.FileType('w'), help="output file name")
parser.add_argument('-i', dest="Col4anchorIDs", type=int, default=10, help="see below")

args = parser.parse_args()


###################################
### 1. read and parse <GO_anno> ###
###################################
col_GOterm = 5
col_GOcategory = 9 # P, F, and C
col_anchorID = args.Col4anchorIDs # default = 10; you may want to use "-i 2" if the <GO_anno> is in .gaf format 

GOterm = 0
GOcategory = "" # P, F, and C
anchorID = ""

GOref_dict = dict() # key = GOcategory (P, F, or C), value = a dict of anchorID and set of GOterms
GOref_dict["P"] = dict() # dict with key = anchorID, value = set of GOterm
GOref_dict["F"] = dict() # dict with key = anchorID, value = set of GOterm
GOref_dict["C"] = dict() # dict with key = anchorID, value = set of GOterm

num_GOref_dict = dict()
num_GOref_dict["P"] = 0
num_GOref_dict["F"] = 0
num_GOref_dict["C"] = 0

num_line = 0

print( "\nreading and parsing %s\nextracting GOterm, GOterm category (P/F/C), and anchorID from columns %d, %d, and %d, respectively," % \
			(args.GO_anno.name, col_GOterm, col_GOcategory, col_anchorID) )

for line in args.GO_anno:
	num_line += 1
	tok = line.split('\t')
	if len(tok) < 10:
		print( "line %d has not enough number of columns, perhaps it's a header line," % num_line )
	else:
		try:
			GOterm = int( tok[col_GOterm - 1].strip().replace("GO:", "") )
			GOcategory = tok[col_GOcategory - 1].strip()
			anchorID = tok[col_anchorID - 1].strip()		
			if anchorID not in GOref_dict[GOcategory]:
				GOref_dict[GOcategory][anchorID] = set() # initializing
			GOref_dict[GOcategory][anchorID].add(GOterm)
			num_GOref_dict[GOcategory] += 1
		except (ValueError, IndexError) as err:
			print( "ignoring line %d in %s : %s" % (num_line, args.GO_anno.name, err) )
			
for key in sorted(GOref_dict, reverse = True):
	print( "for GOterm category %s, %d anchorIDs and %d GOterms captured," % ( key, len(GOref_dict[key]), num_GOref_dict[key] ) )

args.GO_anno.close()


##################################################
### 2. read <anchor_list> and transfer GOterms ###
##################################################

GOnew_dict = dict() # key = GOcategory (P, F, or C), value = a dict of geneID and set of GOterms
GOnew_dict["P"] = dict() # dict with key = geneID, value = set of GOterm
GOnew_dict["F"] = dict() # dict with key = geneID, value = set of GOterm
GOnew_dict["C"] = dict() # dict with key = geneID, value = set of GOterm

num_GOnew_dict = dict()
num_GOnew_dict["P"] = 0
num_GOnew_dict["F"] = 0
num_GOnew_dict["C"] = 0

num_line = 0
geneID = ""
anchorID = ""
spcsName = args.anchor_list.name

print( "\nreading %s and transferring GO annotation," % (args.anchor_list.name) )

# transfer GO annotation
for line in args.anchor_list:
	num_line += 1
	tok = line.split('\t')
	try:
		geneID = tok[0].strip()
		anchorID = tok[1].strip()
		for key in sorted(GOref_dict, reverse = True):
			if anchorID in GOref_dict[key]:
				if geneID not in GOnew_dict[key]:
					GOnew_dict[key][geneID] = set()
				GOnew_dict[key][geneID].update( GOref_dict[key][anchorID] )
	except (ValueError, IndexError) as err:
		print( "ignoring line %d in %s : %s" % (num_line, args.anchor_list, err) )

# count GOterms in GOnew_dict and report 
for key in sorted(GOnew_dict, reverse = True):
	for geneID in GOnew_dict[key]:
		num_GOnew_dict[key] += len( GOnew_dict[key][geneID] )
	
for key in sorted(GOnew_dict, reverse = True):
	print( "for GOterm category %s, %d genes in %s were annotated with %d GOterms," % \
			( key, len(GOnew_dict[key]), args.anchor_list.name, num_GOnew_dict[key] ))

args.anchor_list.close()


#######################
### 3. print output ###
#######################

GOcategory_dict = dict() # key = GOcategory (P, F, or C), value = full name of GOcategory
GOcategory_dict["P"] = "Biological Process"
GOcategory_dict["F"] = "Molecular Function"
GOcategory_dict["C"] = "Cellular Component"

print( "\nwriting output to %s," % args.output.name )

for key in sorted(GOnew_dict, reverse = True):
	# print the header
	args.output.write( "(species=%s)(type=%s)(curator=GO)\n" % ( spcsName, GOcategory_dict[key] ) )
	for geneID in sorted(GOnew_dict[key]):
		for GOterm in sorted(GOnew_dict[key][geneID]):
			args.output.write( "%s = %d\n" % (geneID, GOterm))

args.output.close()

print( "done" )




