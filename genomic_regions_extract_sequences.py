#!/usr/bin/env python
import sys, os, subprocess, argparse
from argparse import RawTextHelpFormatter

###################################################
### 0. script description and parsing arguments ###
###################################################
synopsis1 = "\
  given the .gtfParsed.txt and .genome.fa files, extract intergenic sequences,\n\
  either 5' or 3' regions of all gene models\n"
synopsis2 = "detailed description:\n\
 1. Input files and options:\n\
  - './<Project>.list' including all species IDs (spcsIDs), one per line,\n\
     that appear in OrthNets,\n\
  - '-g Path2Genome': path to genome sequences; expect fasta files named as\n\
     '<spcsID>.genome.fa' in the Path2Genome ['./'],\n\
  - '-t Path2TDfiles': path to 'TDfiles'; 'CL_finder_multi.py -h' for details\n\
  - '-T TDfile_nameFmt': expects TDfiles (or '.gtfParsed.txt' files) named as\n\
     spcsID + TDfile_nameFmt ['.gtfParsed.txt'],\n\
  - '-d': read only IUPAC nucleotide sequences from a ""dirty,"" fasta file,\n\
     i.e. sequence contains spaces, etc. [False].\n\
 2. Extracting and printing intergenic sequences:\n\
  -  the entire intergenic sequences for each genome and lastz commands (-p)\n\
     will be printed on the current folder './'\n\
  - '-l max_len': extracts max_len nucleotides from the 5' (or 3') of each CDS;\n\
     accepts a positive integer; if max_len==0, extracts up to the neighboring\n\
     CDS [0],\n\
  - '-m min_len': do not print extracted sequences shorter than min_len [0],\n\
  - '-r': extracts max_len nucleotides from the 5' (or 3') of each mRNA instead;\n\
     of CDS [False],\n\
  - '-3': extracts 3' intergenic sequences; default is to extract 5' (promoter),\n\
  - '-s': extracted sequences start from the side closer to the gene; e.g. 5'\n\
     sequences start with the start codon or TSS and inverted (default=False; \n\
     i.e. 5' sequences end with the start codon or TSS)\n\
  - '-o Path2temp': path to output and temporary files ['./temp']; CAUTION:\n\
     when run with -p (see below), this folder may contain a large number of files,\n\
 3. Creating files (commands and fasta) for pairwise lastz runs:\n\
  - '-p Pairs2compare': 'Pairs2compare.list' is a tab-delimited file of unique\n\
     pairID (uID), query geneID (qID), and target geneID (tID), a pair per line;\n\
  - with 'Pairs2compare.list' given, print the following:\n\
   (1) a bash file including commands to compare intergenic sequences for all\n\
     pairs as 'Pairs2compare.5p.lastz.sh' (3p instead of 5p with '-3'); this\n\
     will be printed to the current folder\n\
   (2) sequences for all qID and tID are printed as .fa files to Path2temp,\n\
  - '-O options_lastz': optional string passed to lastz; the following is given\n\
     by default: '--chain --seed=111101110010111 --ambiguous=iupac';\n\
     use '-O' to add more options: e.g. -O ' --hspthresh=1500'\n\
  - '-B': use BLASTN output format for lastz (i.e. --format=BLASTN)\n\
  - '-M': use MAF output format for lastz (i.e. --format=maf-)\n\
     * without '-B' or '-M', default output option is 'general', i.e.:\n\
     --format=general:name1,start1,end1,length1,size1,start2,end2,length2,size2\n\
              ,identity,continuity,coverage\n\
  - '-L': print lastz commands only, without extracting sequences [False],\n\
 3.1. Printing fasta files for MSA (multiple sequence alignment):\n\
  - '-P OGs2compare': 'OGs2compare.list' is a tab-delimited file of unique\n\
     pairID (uID), followed by any number of geneIDs; the 5' (or 3') sequences\n\
     for all geneIDs for each uID are printed as 'prefix_uID.fa' in Path2temp,\n\
     where 'prefix_' indicates the regions extracted (e.g. '5p1k_' for 5' 1Kb),\n\
  - '-P OGs2compare' and '-p Pairs2compare' options are mutually exclusive,\n\
 4. Misc:\n\
  - qIDs and tIDs in 'Pairs2compare' (or 'OGs2compare') should be formatted as \n\
     'spcsID|geneID'\n\
 by ohdongha@gmail.com ver0.3.1 20200928\n"
 
#version_history
#20201114 script renamed to "genomic_regions_extract_sequences.py"
#20200928 ver 0.3.1 # activate the '-3' option perhaps by flipping the strand sign? ... hehehe (evil grin) 
#20200614 ver 0.3 # modified to work with python 3.8
#20200311 ver 0.2.1 # minor bug fix for -m option (not printing sequences SHORTER than min_len)
#20200224 ver 0.2 # -m option added
#20190917 ver 0.1.4 # minot bug fix
#20190804 ver 0.1.3 # -M option added to MAF output for lastz; with -p, if 'Pairs2compare' comes in the form of *.list, take the * part as 'Pairs2compare'; minor improvements on lastz commandline output
#20190506 ver 0.1.2 # minor bug fix (fixed prefix when extracting <1kb; add the prefix to lastz output filename); '--format==BLASTN' is NOT default (use -B)  
#20190331 ver 0.1.1 # '--format==BLASTN' will be the default, a counter added ... ; warnings for missing scaffolds 
#20190320 ver 0.1 # use TDfiles instead of CLfm output to get gene_coords; print bash commands for lastz
#20190217 ver 0.0

parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)

# positional arguments
parser.add_argument('Project', type=str, help="'./<Project>.list' includes spcsIDs being compared")
# options and parameters
parser.add_argument('-g', dest="Path2Genome", type=str, default="./", help='["./"]; see below')
parser.add_argument('-t', dest="Path2TDfiles", type=str, default="./", help='["./"]; see below')
parser.add_argument('-T', dest="TDfile_nameFmt", type=str, default=".gtfParsed.txt", help='[".gtfParsed.txt"]; see below')
parser.add_argument('-o', dest="Path2temp", type=str, default="./temp", help="see below")
parser.add_argument('-l', dest="max_len", type=int, default=0, help="see below")
parser.add_argument('-m', dest="min_len", type=int, default=0, help="see below")
parser.add_argument('-d', dest="dirty_seq", action="store_true", default=False)
parser.add_argument('-r', dest="use_mRNA", action="store_true", default=False)
parser.add_argument('-3', dest="three_prime", action="store_true", default=False)
parser.add_argument('-s', dest="start_ATG", action="store_true", default=False)
parser.add_argument('-p', dest="Pairs2compare", type=str, default="__NA__", help="see below")
parser.add_argument('-O', dest="options_lastz", type=str, default="", help="see below")
parser.add_argument('-M', dest="MAF", action="store_true", default=False)
parser.add_argument('-B', dest="BLASTN", action="store_true", default=False)
parser.add_argument('-L', dest="lastz_cmd_only", action="store_true", default=False)
parser.add_argument('-P', dest="OGs2compare", type=str, default="__NA__", help="see below")

args = parser.parse_args()


################################################
### 1. defining global parameters/arguements ###
################################################
# lastz format (when '-p' arguement is given) (modify as needed)
lastz_seed = "111101110010111"

if args.BLASTN:
	lastz_format = "BLASTN"
elif args.MAF:
	lastz_format = "maf-"
else:
	lastz_format = "general:name1,start1,end1,length1,size1,start2,end2,length2,size2,identity,continuity,coverage" #name2 is not needed since pairID is used for the name
lastz_cmd_backbone = "lastz --chain --seed=%s --ambiguous=iupac --format=%s %s " % (lastz_seed, lastz_format, args.options_lastz)

# defining PATHs and create Output directory, if not already exisiting
path_genome = args.Path2Genome
if path_genome[-1] != "/": path_genome = path_genome + "/"
path_TDfiles = args.Path2TDfiles
if path_TDfiles[-1] != "/": path_TDfiles = path_TDfiles + "/"

#define sequence ID prefix (sID_prefix)
sID_prefix = "5p"
if args.three_prime: # define sID_prefix
	sID_prefix = "3p"
if args.max_len >= 1000:
	sID_prefix += ("%.1fk" % (args.max_len/1000.0)).replace(".0", "") # 5p1k, 5p1.5k, 5p2k, etc ...  
elif args.max_len > 0: # 1~999
	sID_prefix += "%d" % args.max_len

# decide whether to print lastz commands, etc, fasta files for MSA, or nothing 
print_lastz_commands = False
print_OG_fasta = False
path_Pairs2compare = ""
path_OGs2compare = ""
if args.Pairs2compare != "__NA__":
	if os.path.isfile( args.Pairs2compare + '.list' ):
		print_lastz_commands = True
		path_Pairs2compare = args.Pairs2compare + '.list'
	elif os.path.isfile( args.Pairs2compare ):
		print_lastz_commands = True
		path_Pairs2compare = args.Pairs2compare
	else:
		print( "to create files for lastz, the following file is needed: %s" % ( args.Pairs2compare + '.list' ) )
		print( "will proceed extracting intergenic sequences, but files for lastz won't be created" )
elif args.OGs2compare != "__NA__":
	if os.path.isfile( args.OGs2compare + '.list' ):
		print_OG_fasta = True
		path_OGs2compare = args.OGs2compare + '.list'
	elif os.path.isfile( args.OGs2compare ):
		print_OG_fasta = True
		path_OGs2compare = args.OGs2compare
	else:
		print( "to create fasta files for multiple sequence alignment (MSA), the following list file is needed: %s" % ( args.OGs2compare + '.list' ) )
		print( "will proceed extracting intergenic sequences, but fasta files for MSA won't be created" )
else:
	print( "Neither '-p' or '-P' arguement is given - will print sequences for all species only," )

# create "temp" folder to output subsets of sequences with '-p' or '-P'	
if print_lastz_commands and print_OG_fasta:
	print( "Error: -P and -p arguments are mutually exclusive; see 'extract_intergenic.py -h'" )
	sys.exit(1)
elif print_lastz_commands or print_OG_fasta: # path_temp required only when either '-P' or '-p' is used.
	path_temp = args.Path2temp
	if path_temp[-1] != "/": path_temp = path_temp + "/"
	try: 
		os.makedirs(path_temp)
	except OSError:
		if not os.path.isdir(path_temp): raise

#function to get nucleotides from a string # use this only if ".genome.fa" is dirty (slow)
def get_nucleotide(str1):
	return "".join(re.findall('[ATGCNatgcnRYSWKMBDHVryswkmbdhv]',str1))	
		
#function to find complementing nucleotides
def invert(seq):
	complement = []
	for c in seq:
		if(c == 'a'): complement.append( 't' )
		elif(c == 'A'): complement.append( 'T' )
		elif(c == 't'): complement.append( 'a' )
		elif(c == 'T'): complement.append( 'A' )
		elif(c == 'g'): complement.append( 'c' )
		elif(c == 'G'): complement.append( 'C' )
		elif(c == 'c'): complement.append( 'g' )
		elif(c == 'C'): complement.append( 'G' )
		elif(c == 'R'): complement.append( 'Y' )
		elif(c == 'r'): complement.append( 'y' )
		elif(c == 'Y'): complement.append( 'R' )
		elif(c == 'y'): complement.append( 'r' )
		elif(c == 'S'): complement.append( 'W' )
		elif(c == 's'): complement.append( 'w' )
		elif(c == 'W'): complement.append( 'S' )
		elif(c == 'w'): complement.append( 's' )
		elif(c == 'K'): complement.append( 'M' )
		elif(c == 'k'): complement.append( 'm' )
		elif(c == 'M'): complement.append( 'K' )
		elif(c == 'm'): complement.append( 'k' )
		elif(c == 'B'): complement.append( 'V' )
		elif(c == 'V'): complement.append( 'B' )
		elif(c == 'b'): complement.append( 'v' )
		elif(c == 'v'): complement.append( 'b' )
		elif(c == 'D'): complement.append( 'H' )
		elif(c == 'd'): complement.append( 'h' )
		elif(c == 'H'): complement.append( 'D' )
		elif(c == 'h'): complement.append( 'd' )
		# Put more rules here!!!
		elif(c != '\n'): complement.append( str(c) )
	complement.reverse()
	return ''.join( complement )

	
################################
### 2. reading the list file ###
################################
projectID = args.Project
try:
	fin_SpcsList = open(projectID + '.list', 'r')
except IOError:
	fin_SpcsList = open(projectID, 'r')
	projectID = projectID[:-5]
spcsID_list = []

print( "\nreading the list file:" + fin_SpcsList.name )
for line in fin_SpcsList:
	spcsID_list.append(line.strip())
	print( line.strip() )
fin_SpcsList.close()
print( "Total %d species IDs detected." % (len(spcsID_list)) )


########################################################################
### 2+. reading the pair file (with '-p' ) and print lastz commands, ###
###  or reading the OG file (with '-P' ) and store ortholog sIDs     ###
########################################################################
if print_lastz_commands:
	fin_pairs = open( path_Pairs2compare, 'r')
	pairs_base = os.path.basename( path_Pairs2compare ) 
	if pairs_base.endswith(".list"):
		pairs_base = pairs_base[:-5] # remove ".list"
	fout_bash = open( pairs_base + ".%s.lastz.sh" % sID_prefix, 'w' ) # print in the current folder

	lastz_output_path = pairs_base + ".%s.lastz.output.txt" % sID_prefix # lastz output and log are printed to the current folder
	lastz_log_path = pairs_base + ".%s.lastz.log" % sID_prefix
		
	print( "\ncreating lastz commands to compare pairs of sequences in %s" % fin_pairs.name )
	print( "expecting in %s unique pairID, query geneID, and target geneID; each geneID is formatted as 'spcsID|geneID'" % fin_pairs.name )
	
	num_line = 0
	pairs_dict = dict() # key = uID, value = [qID, tID]
	sID_2bCompared_set = set() # set of sequence IDs to be compared by lastz
	uID = ""
	qID = ""
	tID = ""
	
	for line in fin_pairs:
		tok = line.split('\t')
		num_line += 1
		skip = False
	
		if len(tok) < 3:
			print( "not enough columns in line %d" % num_line )
		else:
			uID = tok[0].strip()
			qID = tok[1].strip()
			tID = tok[2].strip()
			if len( qID.split('|') ) != 2 or len( tID.split('|') ) != 2 :
				print( "illegal pair of geneIDs line %d (perhaps it's a header?), skipping: %s" % (num_line, line.strip()) )
			else:
				pairs_dict[ uID ] = [ qID, tID ]
				sID_2bCompared_set.add( sID_prefix + '_' + qID )
				sID_2bCompared_set.add( sID_prefix + '_' + tID )
				
				query_path = ( path_temp + "%s.fa" % (sID_prefix + '_' + qID) ).replace('|', ':') # having '|' in a file name is not a good idea -_-;;;
				target_path = ( path_temp + "%s.fa" % (sID_prefix + '_' + tID) ).replace('|', ':')

				query_path = uID + "_seq1::" + query_path # asking lastz to use the uID as "nicknames" for both query and target
				target_path = uID + "_seq2::" + target_path
				
				if num_line == 1:	
					lastz_cmd = lastz_cmd_backbone + "%s %s > %s 2> %s\n" % ( query_path, target_path, lastz_output_path, lastz_log_path) # create
				else:
					lastz_cmd = lastz_cmd_backbone + "%s %s >> %s 2>> %s\n" % ( query_path, target_path, lastz_output_path, lastz_log_path) # append
				fout_bash.write( lastz_cmd )
	fin_pairs.close()
	
elif print_OG_fasta:
	fin_OGs = open( path_OGs2compare, 'r')
	orthologs_dict = dict() # key = uID, value = list of ortholog geneIDs
	sID_2bCompared_set = set() # set of sequence IDs to be compared by lastz

	uID = ""
	orthologID_list = []

	fin_OGs = open( path_OGs2compare, 'r')
	
	num_line = 0
	orthologs_dict = dict() # key = uID, value = list of ortholog geneIDs

	uID = ""
	orthologID_list = []
	
	for line in fin_OGs:
		tok = line.split('\t')
		num_line += 1
		skip = False
	
		if len(tok) < 3:
			print( "not enough columns in line %d, skipping; at least 2 sequence IDs needed for MSA," % num_line )
		else:
			uID = tok[0].strip()
			orthologID_list = tok[1:]
			sID_check = True
			for s in orthologID_list:
				if len( s.strip().split('|') ) != 2 :
					print( "illegal sequence ID: '%s' found in line %d (header?), skipping: %s" % (s.strip(), num_line, line.strip()) )
					sID_check = False
					break
				else:
					sID_2bCompared_set.add( sID_prefix + '_' + s.strip() )
			
			if sID_check:
				orthologs_dict[ uID ] = orthologID_list
			
	fin_OGs.close()


##########################################
### 3. extracting intergenic sequences ###
##########################################
## 3.0 define global dict, arguments, and parameters
if not args.lastz_cmd_only:
	chr_seq_dict = dict() # key = chrID (cID); value = sequence as a string
	chr_len_dict = dict() # key = chrID (cID); value = length of the chr/scf/contig
	gene_coords_dict = dict() # key = CLfm num_line, value = [geneID (gID), cID, Str, CDS/mRNA_s, CDS/mRNA_e]
	coords_2cut_dict = dict() # key = CLfm num_line, value = [seqID (sID), cID, Str, start, end]
	
	seq_2bCompared_dict = dict() # key = sID, value = sequence # this will be used in ## 2.4 with '-p' option
	seqHeader_dict = dict() # key = sID, value = sequence header line, including the Chr coordinates, etc
	
	if args.use_mRNA: # which columns to use in the CLfm files
		index_2use = [0, 1, 2, 3, 4]
	else:
		index_2use = [0, 1, 2, 7, 8]
	
		
	## start iterating over spcsID_list
	for spcsID in spcsID_list:
		## 3.1 reading .genome.fa file
		fin_genome = open( path_genome + spcsID + ".genome.fa", 'r')
		print( "\nreading genome sequences from %s" % fin_genome.name )
		
		chr_seq_dict.clear() # initialize
		cID=""
		seq_in_line = ""
		collected_seq = ""	
		seq_count = 0
		nt_count = 0
		collect = False
		
		for line in fin_genome:
			if line[0] == '>':
				if cID != "" : # for the first sequence
					chr_seq_dict[cID] = collected_seq
				cID = line[1:-1].split()[0]
	#			print cID # for debugging purpose
				collected_seq = ""
				seq_in_line = ""
				seq_count = seq_count + 1
				collect = True
			elif collect :
				if args.dirty_seq:
					seq_in_line = get_nucleotide(line)  ## slow; only if you expect dirty fasta
				else:
					seq_in_line = line.strip()
				nt_count = nt_count + len(seq_in_line)
				collected_seq = collected_seq + seq_in_line
				
		chr_seq_dict[cID] = collected_seq # for the last sequence
		print( "total %d sequences, %d nucleotides read from %s ... " % \
				(seq_count, nt_count, fin_genome.name) )
		fin_genome.close()
	
		chr_len_dict.clear() # seq_len dictionary
		for cID in chr_seq_dict:
			chr_len_dict[cID] = len( chr_seq_dict[cID] )
		
	
		## 3.2 reading gene coordinates from the CLfm file
		TDfileName = path_TDfiles + spcsID + args.TDfile_nameFmt
		fin_TDfile = open( TDfileName, 'r')
		print( "reading gene model coordinates from %s" % fin_TDfile.name )
	
		gene_coords_dict.clear() # initialize
		n = 0
		header = True
		
		for line in fin_TDfile:
			tok = line.strip().split('\t')
			if header:
				header = False
			else:
				n += 1
				gene_coords_dict[n] = [ tok[i] for i in index_2use ]
		print( "total %d gene model coordinates read from %s" % ( n, fin_TDfile.name) )
		fin_TDfile.close()
		
	
		## 3.3 obtain coordinates to extract
		#gene_coords_dict # key = CLfm num_line, value = [gID, cID, Str, CDS/mRNA_s, CDS/mRNA_e]
		#coords_2cut_dict # key = CLfm num_line, value = [sID, cID, Str, start, end]	
		coords_2cut_dict.clear() # initialize
		sID = ""
		cID = ""
		strand = ""
		gene_s = 0
		gene_e = 0
		coord_s = 0
		coord_e = 0
		
		for n in sorted( gene_coords_dict ): # iterate over genes as they appear in the CLfm file
			sID = sID_prefix + '_' + spcsID + '|' + gene_coords_dict[n][0] # example sID = "5p1k_spcsID|geneID"
			cID = gene_coords_dict[n][1]
			strand = gene_coords_dict[n][2]
			gene_s = int( gene_coords_dict[n][3] )
			gene_e = int( gene_coords_dict[n][4] )
		
			coords_2cut_dict[n] = [sID, cID, strand] # first three records
			
			if args.three_prime: # let's try just reversing strand?
				if strand == '+':
					strand = '-'
				elif strand == '-':
					strand = '+' 		
			if cID in chr_len_dict:
				if strand == '+':
					coord_e = min( gene_s + 2, chr_len_dict[ cID ] ) # include "ATG" 
					if args.max_len == 0 : # extract up to the end of the previous gene model
						coord_s = 1
						if ( n-1 ) in gene_coords_dict:
							if gene_coords_dict[n-1][1] == cID: 
								coord_s = int( gene_coords_dict[n-1][4] ) + 1
					else:
						coord_s = max( gene_s - args.max_len, 1 )
				elif strand == '-':
					coord_s = max( gene_e - 2, 1 ) # include "ATG"
					if args.max_len == 0 : # extract up to the start of the next gene model
						coord_e = chr_len_dict[ cID ]
						if ( n+1 ) in gene_coords_dict:
							if gene_coords_dict[n+1][1] == cID: 
								coord_e = int( gene_coords_dict[n+1][3] ) - 1
					else:
						coord_e = min( gene_e + args.max_len, chr_len_dict[ cID ])
				else:
					print( "Warning: ambiguous strand chracter %s for seqID: %s, skipping" % ( strand, sID ) )   
			else:
				print( "Warning: scaffold/chromosome: %s is not in the genome, skipping geneID: %s" % ( cID, sID ) )  
		
			coords_2cut_dict[n] = coords_2cut_dict[n] + [ coord_s, coord_e ]
		
		## 3.4 extract and print sequences (with '-p' or '-P' also store sequences in a dictionary) 
		out_fileName = "./" + sID_prefix + '_' + spcsID + ".fa" # print to the current folder
		fout = open( out_fileName, 'w')
		print( "writing extracted sequences to %s" % fout.name )
	
		sID = ""
		cID = ""
		strand = ""
		coord_s = 0
		coord_e = 0
		seq_extracted = ""
			
		for n in sorted( coords_2cut_dict ): # iterate over genes as they appear in the CLfm file
			sID = coords_2cut_dict[n][0] # example sID = "5p1k_spcsID|geneID"
			cID = coords_2cut_dict[n][1]
			strand = coords_2cut_dict[n][2]
			coord_s = int( coords_2cut_dict[n][3] )
			coord_e = int( coords_2cut_dict[n][4] )
			flip = False
			header_line = ""
	
			if cID in chr_seq_dict:
				if coord_s >= 0 and coord_e <= chr_len_dict[cID]:
					seq_extracted = chr_seq_dict[cID][coord_s-1 : coord_e]
					header_line = ">%s\t%s\t%s" % ( sID, \
								'|'.join( gene_coords_dict[n] ), \
								cID + ':' + str(coord_s) + '-' + str(coord_e))
				if args.three_prime: # decide whether to flip
					if ( args.start_ATG and strand == '-' ) or ( (not args.start_ATG) and strand == '+' ):
						flip = True
				else:
					if ( args.start_ATG and strand == '+' ) or ( (not args.start_ATG) and strand == '-' ):
						flip = True
						
				if flip:
					seq_extracted = invert( seq_extracted )
					header_line += "_inv"
				
				if len(seq_extracted) >= args.min_len: # v0.2 do not print if shorter than min_len
					header_line = header_line + "\tlen=%d" % len(seq_extracted)
					fout.write( header_line + '\n' + seq_extracted + '\n' )
				else:
					print( "skipping too short sequence: %s in %s" % ( header_line, fout.name) )
					
				if ( print_lastz_commands or print_OG_fasta ) and sID in sID_2bCompared_set:
					seq_2bCompared_dict[ sID ] = seq_extracted
					seqHeader_dict[ sID ] = header_line[1:] # let's not include ">"
			else:
				print( "Warning: cID %s for sID %s is not found in the genome.fa, skipping" % ( cID, sID ) )
	
		print( "done extracting intergenic sequences for %s" % spcsID )
		fout.close()
	

##########################################################################################
### 3+. with '-p' (or '-P') print sequences to be used by lastz (or other MSA program) ###
##########################################################################################
if print_lastz_commands and not args.lastz_cmd_only: 		
	print( "\nprinting sequences as individual fasta files ready for lastz to %s" % path_temp )
	print( "CAUTION: %s may contain a large number of .fasta files (one file per each sequence)..." % path_temp )
	seqID = ""
	seqPath = ""
	already_printed = set()
	
	total_pairs = len( pairs_dict )
	num_processed = 0
	
	for uID in sorted( pairs_dict ):
		for i in [0, 1]:
			num_processed += 1
			seqID = sID_prefix + '_' + pairs_dict[ uID ][i]
			if seqID not in already_printed :
				seqPath = ( path_temp + "%s.fa" % (seqID) ).replace('|', ':')
				if not os.path.isfile( seqPath ): # skip if the sequence file already exist
					if seqID in seq_2bCompared_dict:
						if len( seq_2bCompared_dict[ seqID ] ) > args.min_len:
							fout_seq = open( seqPath, 'w')
							fout_seq.write( '>%s\n' % seqID )
							fout_seq.write( seq_2bCompared_dict[ seqID ] + '\n' )
							fout_seq.close()
						else:
							print( "skipping too short sequence: %s\tlen=%d in %s" % ( seqID, len( seq_2bCompared_dict[ seqID ] ), seqPath) )
					else:
						print( "sequence for the ID %s is not available, perhaps it is part of the header line?" % seqID )
				already_printed.add( seqID )
				
		# counter display
		if ( num_processed % 10000 == 0):
			sys.stdout.write("\r   processed %d / %d pairs; %d sequences printed to %s " \
				% (num_processed, total_pairs, len( already_printed ), path_temp) )
			sys.stdout.flush()

	print( "\ndoen printing fasta files for lastz run, \nnow good luck running the bash file: %s" % fout_bash.name )
	fout_bash.close()

elif print_OG_fasta: 		
	print( "\nprinting intergenic sequence fasta files for genes listed in %s for multiple sequence alignment (MSA) to %s" % (path_OGs2compare, path_temp) )
	print( "CAUTION: %s may contain a large number of .fasta files (one file per each ortholog group)..." % path_temp )
	seqID = ""
	seqPath = ""
	
	total_pairs = len( orthologs_dict )
	num_processed = 0

	for uID in sorted( orthologs_dict ):
		num_processed += 1
		seqPath = ( path_temp + "%s.fa" % (sID_prefix + '_' + uID) )
		fout_seq = open( seqPath, 'w')
	
		for s in orthologs_dict[ uID ]:
			seqID = sID_prefix + '_' + s.strip()
			if seqID in seq_2bCompared_dict:
				if len(seq_2bCompared_dict[ seqID ] ) > args.min_len:
					fout_seq.write( '>%s\n' % seqHeader_dict[ seqID ] )
					fout_seq.write( seq_2bCompared_dict[ seqID ] + '\n' )
				else:
					print( "skipping too short sequence: %s in %s" % ( seqHeader_dict[ seqID ], fout_seq.name) )
			else:
				print( "sequence for the ID %s is not available, perhaps it is part of the header line?" % seqID )
		fout_seq.close()
		
		# counter display
		if ( num_processed % 1000 == 0):
			sys.stdout.write("\r   wrote fasta files for %d / total %d ortholog groups to %s " \
				% (num_processed, total_pairs, path_temp) )
			sys.stdout.flush()

	print( "\ndoen printing fasta files for MSA in %s, \nnow good luck running the a multi-sequence aligner, such as t_coffee, muscle, or fsa on them," % path_temp )
		