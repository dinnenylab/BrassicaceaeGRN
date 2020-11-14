#!/usr/bin/env python
import sys, os, argparse
from argparse import RawTextHelpFormatter

##############################
### 0.1 script description ###
##############################
synopsis1 = "\
  - perform multiple test correction on a column of p-values;\n\
  - currently support:\n\
     'fdr_bh' and other methods from statsmodels.stats.multitest\n\
     'qvalue' method from multipy ( https://puolival.github.io/multipy/ )\n"
  
synopsis2 = "detailed description:\n\
 1. Input files:\n\
  - <input>: tab-delimited file with a header and p-values in the Nth column;\n\
     all non-numerical values in the Nth column is considered p = 1.0;\n\
  - <N>: col. number that includes p-values to apply multiple test correction;\n\
     N=1 to use the first column; [1]\n\
  - <output>: output file name; [<input> + '.adjusted']\n\
 2. Options and parameters:\n\
  - '-m method': string for methods from statsmodels.stats.multitest ['fdr_bh']\n\
  - '-q' | '--qvalue': import multipy and print q-values using 'qvalue' method\n\
  - '-e'|'--scientific': adjusted p-values are printed in scientific notations \n\
     (e.g. 1.00e-2) [False];\n\
 3. Output:\n\
  - print the adjusted p-values (or q-values with '\q') as the last column;\n\
  - add '.adj' (or '.q' with '-q') to the header of Nth column and use as the new\n\
     column name;\n\
by ohdongha@gmail.com 2020623 ver 0.0\n\n"
#version_history
#20201111 ver 0.0.1 '-e' option to print adjusted p-values in scientific annotations
#20200623 ver 0.0


#############################
### 0.2 parsing arguments ###
#############################
parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)
# positional arguments
parser.add_argument('input', type=argparse.FileType('r'))
parser.add_argument('N', type=int, default=1)
parser.add_argument('output', type=str, default="__NA__")
# options
parser.add_argument('-m', dest="method", type=str, default="fdr_bh")
parser.add_argument('-q', '--qvalue', action="store_true", default=False) 
parser.add_argument('-e', '--scientific', action="store_true", default=False) 

# parsing arguments
args = parser.parse_args()

try:
	colNumber = int( args.N ) - 1
	if colNumber < 0:
		raise ValueError
except ValueError:
	print("Warning: invalid value for the column number N, using default N=1")
	colNumber = 0

if args.output == "__NA__":
#	outfile_base = os.path.splitext(args.input.name)[0] # add "_pplt.png" etc later
	outfile_base = args.input.name + ".adjusted"
else:
	outfile_base = args.output
	

###########################
### 0.3 importing stuff ###
###########################
try:
	if args.qvalue:
		from multipy.fdr import qvalue
	else:
		import statsmodels.stats.multitest as multi
except ImportError as ErrorMessage:
	print(str(ErrorMessage))
	sys.exit()


######################################################
### 1. read input and correct for multiple testing ###
######################################################
## 1.1 reading input
header = True
colName = ""
lineNumber = 0
Pval_list = [] # list of p-values

for line in args.input:
	tok = line.split('\t')
	lineNumber += 1
	Pval = 1e+0
	
	if header:
		colName = tok[ colNumber ]
		header = False
	else:
		try:
			Pval = float( tok[ colNumber ] )
		except ValueError:
			print("Warning: line #%d does not have a number in column #%d, consider p as 1.0 " % (lineNumber, colNumber + 1) )
			Pval = 1.0
		Pval_list.append( Pval )

## 1.2 correct for multiple testing
Padj_list = [] # list of adjusted p-values (or q-values)

if args.qvalue:
	Padj_list = list( qvalue(Pval_list, verbose = False)[1] )
else:
	Padj_list = list( multi.multipletests(Pval_list, method = args.method)[1] )


#######################
### 2. print output ###
#######################
fout = open( outfile_base, 'w')

args.input.seek(0)
header = True
lineNumber = 0

for line in args.input:
	if header:
		if args.qvalue:
			fout.write( line.strip() + "\t%s.q\n" % colName )
		else:
			fout.write( line.strip() + "\t%s.adj\n" % colName )
		header = False
	else:
		try: 
			if args.scientific:
				fout.write( line.strip() + '\t%.2e\n' % Padj_list[ lineNumber ] )			
			else:
				fout.write( line.strip() + '\t%f\n' % Padj_list[ lineNumber ] )
		except IndexError:
			print("lineNumber = %d, len(Padj_list) = %d" % (lineNumber, len(Padj_list) ) )
		lineNumber += 1

print("## done writing adjust p-values (or q-values) as the last column of %s." % fout.name)

fout.close()
args.input.close()

