#!/usr/bin/env python
import sys, os, math, subprocess, datetime, argparse
from argparse import RawTextHelpFormatter


####################################################
### 0.1 script description and parsing arguments ###
####################################################
synopsis1 = "\
  - read protein-coding gene models from a .gtf file and identify genomic\n\
     coordinates for proximal/distal 5' and 3' intergenic regions, exons,\n\
     introns, etc. (i.e. annotating genomic regions) "
synopsis2 = "detailed description:\n\
 1. Input files:\n\
  - <input_gtf> is a .gtf containing CDS records for transcript (or gene) models;\n\
     exon records are optional (under development)\n\
  - <input_genome_index> is a tab-delimited file with ID (chrID) and size for each\n\
     chromosome, scaffold, or contig, one per line; this will be used to prevent\n\
     the genomic regions overflow on a chromosome, etc.; can be obtained by\n\
        count_nt2.py genome.fa genome\n\
     the output file 'genome.list' can be used as <input_genome_index>,\n\
 2. parsing a .gtf file and annotate genomic regions:\n\
  - from <input.gtf>, find out coordinates of CDS and annotate genomic regions\n\
     as one of the following:\n\
	 (1) 5pD: 5' distal (-2000 ~ -501bp from the start of an ORF),\n\
     (2) 5pP: 5' proximal (-500 ~ -1 bp from the start of an ORF),\n\
	 (3) 3pD: 3' distal (+501 ~ +2000bp from the end of an ORF),\n\
     (4) 3pP: 3' proximal (+1 ~ +500bp from the end of an ORF),\n\
     (5) Exo: Exons of an ORF,\n\
     (6) IntF: The first and last intron of an ORF, and\n\
     (7) Int: Introns, except for the first one, of an ORF.\n\
     (8) (optional) 5' and 3' UTRs (under development)\n\
  - works best with a .gtf containing representative protein-coding gene models\n\
     only, with 'CDS' records; 'exon' records are optional (under devleopment),\n\
  - for every CDS features, print all items listed above even if they overlap,\n\
     e.g. a position in a genome can be in the 5' distal intergenic region of\n\
     one gene and 3' proximal region of the other,\n\
 3. Options and parameters:\n\
  - '-g|--gene_id': use 'gene_id' record instead of 'transcript_id' to indicate\n\
     each gene model; [False]\n\
  - '-d|--distal': boundaries for distal intergenic regions; [2000]\n\
  - '-p|--proximal': boundaries for proximal intergenic regions; [500]\n\
  - '-L|--count_last': count 1st and last exons and introns separately (as ExoF,\n\
     ExoL, IntF, and IntL), if there are more than one for an ORF; 1st and last\n\
     exon/introns are not included in Exo and Int categories; [False]\n\
  - '-U|--UTR': (will be added in the future)\n\
 4. Output:\n\
  - <outfile> is tab-delimited, with following fields, one region per line:\n\
     unique ID (uID), chrID, start, end, and depth\n\
  -	 'uID' includes the category (see above), genomic region coordinates, ID of\n\
     the gene to which the genomic region belongs to, and the coordinates of\n\
     the gene,\n\
  -  'depth' can be used to priotize one category over another, when a genomic\n\
     position is annotated with multiple different categories,\n\
by ohdongha@gmail.com 20191011 ver 0.1.1\n\n"

#version_history
#20201114 script renamed to "genomic_regions_annotate.py"
#20191011 ver 0.1.1 when a region is entirely out of bound of a chromosome, print "out_of_bound" as the "depth", 
#20190825 ver 0.1 added '-L' option 
#20190819 ver 0.0 modified from parse_gtf_2table.py 

parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)

# positional parameters
parser.add_argument('input_gtf', type=argparse.FileType('r'))
parser.add_argument('input_genome_index', type=argparse.FileType('r'))
parser.add_argument('outfile', type=argparse.FileType('w'))

# options
parser.add_argument('-g', '--gene_id', action="store_true", default=False)
parser.add_argument('-d', dest="distal", type=int, default= 2000) 
parser.add_argument('-p', dest="proximal", type=int, default= 500) 
parser.add_argument('-U', '--UTR', action="store_true", default=False) # under development
parser.add_argument('-L', '--count_last', action="store_true", default=False)

args = parser.parse_args()
outfile_name = args.outfile.name


#########################################
### 0.2 defining categories and depth ###
#########################################
# defining the depth of categories for each category, feel free to modify
depth_dict = dict() # key = category, value = depth, the higher the more relevant 
depth_dict["5pP"] = 9 # 5pP: 5' proximal (-500 ~ -1 bp from the start of an ORF),
depth_dict["5pD"] = 8 # 5pD: 5' distal (-2000 ~ -501bp from the start of an ORF),
depth_dict["IntF"] = 7 # IntF: The first intron of an ORF,
depth_dict["ExoF"] = 6 # ExoF: The first exons of an ORF, when there more than one ('-L' option) 
depth_dict["3pP"] = 5 # 3pD: 3' distal (+501 ~ +2000bp from the end of an ORF),
depth_dict["3pD"] = 4 # 3pP: 3' proximal (+1 ~ +500bp from the end of an ORF),
depth_dict["IntL"] = 3 # IntL: The last intron of an ORF, when there more than one ('-L' option),
depth_dict["ExoL"] = 2 # ExoL: The last exon of an ORF, when there more than one ('-L' option) 
depth_dict["Int"] = 1 # Int: Introns, except for the first (and the last) one(s), of an ORF.
depth_dict["Exo"] = 0 # Exo: Exons(, except for the first and the last ones) of an ORF, () applies only with '-L' option

# defining dictionaries for all categories of genomic regions
region_dict = dict() # key = ["5pP", "5pD", "IntF", "3pP", "3pD", "Int", "Exo"], value = (see below)

region_dict["Exo"] = dict() # key = geneID; value = dict with key = start and value = end coordinates, for all "exon" records,
region_dict["Int"] = dict() # key = geneID; value = dict with key = start and value = end coordinates, for all "CDS" records,
region_dict["IntF"] = dict() # key = geneID; value = dict with key = start and value = end coordinates, for the 1st intron, when there are more than one,

if args.count_last:
	region_dict["IntL"] = dict() # key = geneID; value = dict with key = start and value = end coordinates, for the 1st intron, when there are more than one,
	region_dict["ExoF"] = dict() # key = geneID; value = dict with key = start and value = end coordinates, for the 1st exon, when there are more than one,
	region_dict["ExoL"] = dict() # key = geneID; value = dict with key = start and value = end coordinates, for the 1st exon, when there are more than one,
	
region_dict["5pD"] = dict() # key = geneID; value = [start, end]
region_dict["3pD"] = dict()
region_dict["5pP"] = dict()
region_dict["3pP"] = dict()


##################################
### 1.1 reading in <input.gtf> ###
##################################
chr_dict = dict() # key = geneID, value = chrID
str_dict = dict() # key = geneID, value = direction of the gene model (+ or -)

## these are reserved for later
#mRNA_start_dict = dict() # key = geneID; value = start of the gene model in the chrID based on "exon" records, the smaller of the two boundaries,
#mRNA_end_dict = dict() # key = geneID; value = end of the gene model in the chrID based on "exon" records, the larger of the two boundaries,

CDS_start_dict = dict() # key = geneID; value = start of the gene model in the chrID based on "CDS" records, the smaller of the two boundaries,
CDS_end_dict = dict() # key = geneID; value = end of the gene model in the chrID based on "CDS" records, the larger of the two boundaries,

newline_accepted = False
chr = ""
type = ""
start = 0
end = 0 
strand = ""
ninthColumn_records = ""
geneID = ""

print "reading %s as the <input.gtf>:" % args.input_gtf.name

for line in args.input_gtf:
	newline_accepted = False
	tok = line.replace('\"','').split('\t')
	try:
		chr = tok[0]
		type = tok[2]
		start = int(tok[3])
		end = int(tok[4])
		strand = tok[6]
		ninthColumn_records = tok[8].split(';')
		for record in ninthColumn_records:
			if args.gene_id:
				if record.strip().split(' ')[0] == 'gene_id':
					geneID = record.strip().split(' ')[1]
					newline_accepted = True
			else:
				if record.strip().split(' ')[0] == 'transcript_id':
					geneID = record.strip().split(' ')[1]
					newline_accepted = True
	except (ValueError, IndexError) :
		print "an invalid line: %s" % line
		newline_accepted = False
	if newline_accepted and type == "CDS" :
		if geneID not in CDS_start_dict:
			chr_dict[geneID] = chr
			str_dict[geneID] = strand
			region_dict["Exo"][geneID] = dict() # initializing
			CDS_start_dict[geneID] = start
			CDS_end_dict[geneID] = end
		else:
			CDS_start_dict[geneID] = min(start, CDS_start_dict[geneID]) # updating the gene boundary 
			CDS_end_dict[geneID] = max(end, CDS_end_dict[geneID])
		region_dict["Exo"][geneID][start] = end # adding a CDS exon
#	elif newline_accepted and type == "exon" : # visit later 
#		if geneID not in mRNA_start_dict:
#			chr_dict[geneID] = chr
#			str_dict[geneID] = strand
#			mRNA_start_dict[geneID] = start
#			mRNA_end_dict[geneID] = end
#			...
#		else:
#			mRNA_start_dict[geneID] = min(start, mRNA_start_dict[geneID])
#			mRNA_end_dict[geneID] = max(end, mRNA_end_dict[geneID])
#			...

print "## %d gene models were found in %s.\n" % ( len(CDS_start_dict) , args.input_gtf.name )
args.input_gtf.close()


###########################################
### 1.2 reading in <input_genome_index> ###
###########################################
chr_len_dict = dict() # key = chrID, value = chr_len

print "reading %s as the <input_genome_index>:" % args.input_genome_index.name
for line in args.input_genome_index:
	tok = line.split('\t')
	try:
		chr_len_dict[ tok[0].split()[0].strip() ] = int( tok[1] ) # only the string before the first space considered as a chrID, 
	except (ValueError, IndexError):
		print "Invalid line found in %s: %s" % (args.input_genome_index.name, line)

print "## %d chromosome/scaffold/contigs read from %s.\n" % ( len( chr_len_dict ), args.input_genome_index.name)
args.input_genome_index.close()


#####################################
### 2. annotating genomic regions ###
#####################################
for g in CDS_start_dict: # g == geneID
	# if there are more than one exons, records introns
	if len( region_dict["Exo"][g] ) >= 2: 
		region_dict["Int"][g] = dict() # key = start and value = end coordinates, for each intron
		exon_end_prev = 0
		for s in sorted( region_dict["Exo"][g] ): # annotating introns
			if exon_end_prev == 0: # i.e. first exon
				exon_end_prev = region_dict["Exo"][g][s]
			else:
				region_dict["Int"][g][ exon_end_prev + 1 ] = s - 1 # intron = (exon_end_prev +1) ~ ( exon_start - 1 ) 
				exon_end_prev = region_dict["Exo"][g][s] # update exon_end_prev
				
	# finding 5' and 3' intergenic regions, as well as "1st intron" when available
	chr_end = chr_len_dict[ chr_dict[g] ]	
	intF_start = -1 # initialize
	if args.count_last: # added at v0.1
		exoF_start = -1 # initialize when '-L'
		exoL_start = -1 
		intL_start = -1 
	if str_dict[g] == '+': 
		region_dict["5pD"][g] = [ max(1, CDS_start_dict[g] - args.distal ), max(1, CDS_start_dict[g] - args.proximal - 1 ) ]
		region_dict["3pD"][g] = [ min(chr_end, CDS_end_dict[g] + args.proximal + 1 ), min(chr_end, CDS_end_dict[g] + args.distal)]
		region_dict["5pP"][g] = [ max(1, CDS_start_dict[g] - args.proximal ), max(1, CDS_start_dict[g] - 1 ) ]
		region_dict["3pP"][g] = [ min(chr_end, CDS_end_dict[g] + 1 ), min(chr_end, CDS_end_dict[g] + args.proximal)]
		if len( region_dict["Exo"][g] ) >= 3: # i.e. if there are more than one introns
			intF_start = min( region_dict["Int"][g] )
			if args.count_last:
				intL_start = max( region_dict["Int"][g] )
		if args.count_last and len( region_dict["Exo"][g] ) >= 2: # i.e. if there are more than one exons
			exoF_start = min( region_dict["Exo"][g] )
			exoL_start = max( region_dict["Exo"][g] )
	elif str_dict[g] == '-':
		region_dict["3pD"][g] = [ max(1, CDS_start_dict[g] - args.distal ), max(1, CDS_start_dict[g] - args.proximal - 1 ) ]
		region_dict["5pD"][g] = [ min(chr_end, CDS_end_dict[g] + args.proximal + 1 ), min(chr_end, CDS_end_dict[g] + args.distal)]
		region_dict["3pP"][g] = [ max(1, CDS_start_dict[g] - args.proximal ), max(1, CDS_start_dict[g] - 1 ) ]
		region_dict["5pP"][g] = [ min(chr_end, CDS_end_dict[g] + 1 ), min(chr_end, CDS_end_dict[g] + args.proximal)]
		if len( region_dict["Exo"][g] ) >= 3: # i.e. if there are more than one introns
			intF_start = max( region_dict["Int"][g] )
			if args.count_last:
				intL_start = min( region_dict["Int"][g] )
		if args.count_last and len( region_dict["Exo"][g] ) >= 2: # i.e. if there are more than one exons
			exoF_start = max( region_dict["Exo"][g] )
			exoL_start = min( region_dict["Exo"][g] )
	# record 1st intron etc
	if intF_start != -1:
		region_dict["IntF"][g] = dict()
		region_dict["IntF"][g][intF_start] = region_dict["Int"][g][intF_start] # use the key intF_start to separate 1st introns from the remaining later ...
		if args.count_last:
			region_dict["IntL"][g] = dict()
			region_dict["IntL"][g][intL_start] = region_dict["Int"][g][intL_start]
	if args.count_last and exoF_start != -1:
		region_dict["ExoF"][g] = dict()
		region_dict["ExoF"][g][exoF_start] = region_dict["Exo"][g][exoF_start]
		region_dict["ExoL"][g] = dict()
		region_dict["ExoL"][g][exoL_start] = region_dict["Exo"][g][exoL_start]


############################
### 3. writing <outfile> ###
############################
categoty = ""
uID = ""
start = 0
end = 0

print "writing annotated genomic regions to %s:" % outfile_name
#args.outfile.write("uID\tchrID\tstart\tend\tannotation\tdepth\n") # writing the header
args.outfile.write("uID\tchrID\tstart\tend\tdepth\n") # writing the header

# printing 5' and 3' intergenic regions 
for c in ["5pP", "5pD", "3pP", "3pD"]: # c == category
	for g in sorted(region_dict[c]): # g == geneID
		start = region_dict[c][g][0]
		end = region_dict[c][g][1]

		uID = "%s|%s:%d-%d_of_%s@%s:%d-%d(%s)" % (c, chr_dict[g], start, end, \
			g, chr_dict[g], CDS_start_dict[g], CDS_end_dict[g], str_dict[g])
			
		# v0.1.1 process 5pD, 3pD, 5pP, and 3pP that are completely out of boundary
		if start == end:
			args.outfile.write( "%s\t%s\t%d\t%d\tout_of_bound\n" % (uID, chr_dict[g], start, end ) )		
		else:
			args.outfile.write( "%s\t%s\t%d\t%d\t%d\n" % (uID, chr_dict[g], start, end, depth_dict[c]) )		
#		uID = "%s_%s|%s:%d-%d" % (c, g, chr_dict[g], start, end)  
#		annotation = "%s_of_%s@%s:%d-%d(%s)" % (c, g, chr_dict[g], CDS_start_dict[g], CDS_end_dict[g], str_dict[g])
#		args.outfile.write( "%s\t%s\t%d\t%d\t%s\t%d\n" % (uID, chr_dict[g], start, end, annotation, depth_dict[c]) )		

# printing exons and introns
skip = False
exon_categories = ["Exo", "IntF", "Int"]
if args.count_last:
	exon_categories = exon_categories + ["IntL", "ExoF", "ExoL"]
for c in exon_categories: # c == category
	for g in sorted(region_dict[c]): # g == geneID
		for s in region_dict[c][g]:
			skip = False

			if c == "Int" and g in region_dict["IntF"]:
				if s in region_dict["IntF"][g]: # among introns, do not print 1st introns twice
					skip = True
			if not skip:
				start = s
				end = region_dict[c][g][s]
	
				uID = "%s|%s:%d-%d_of_%s@%s:%d-%d(%s)" % (c, chr_dict[g], start, end, \
					g, chr_dict[g], CDS_start_dict[g], CDS_end_dict[g], str_dict[g])
				args.outfile.write( "%s\t%s\t%d\t%d\t%d\n" % (uID, chr_dict[g], start, end, depth_dict[c]) )

#				uID = "%s_%s|%s:%d-%d" % (c, g, chr_dict[g], start, end)  
#				annotation = "%s_of_%s@%s:%d-%d(%s)" % (c, g, chr_dict[g], CDS_start_dict[g], CDS_end_dict[g], str_dict[g])
#				args.outfile.write( "%s\t%s\t%d\t%d\t%s\t%d\n" % (uID, chr_dict[g], start, end, annotation, depth_dict[c]) )		

print "## genomic regions for %d protein-coding gene models (ORFs) written to %s.\n" % ( len( CDS_start_dict ), args.outfile.name)
args.outfile.close()

## sort the output file
print "sorting %s:" % outfile_name
i = datetime.datetime.now()
temp_filename = "genomic_regions_" + i.strftime('%Y%m%d_%H%M%S')
subprocess.call("awk 'NR == 1; NR > 1 {print $0 | \"sort -k2,2 -k3,3n\"}' " + outfile_name + " > " + temp_filename , shell=True)
subprocess.call("mv " + temp_filename + " " + outfile_name , shell=True)

print "all done\n"