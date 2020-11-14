#!/usr/bin/env python
import sys, os, argparse
from argparse import RawTextHelpFormatter

##############################
### 0.1 script description ###
##############################
synopsis1 = "\
  - read a table with values in 2x2 contingency tables, one table per row;\n\
  - add chi^2 (or g), p-value, dof, and expected array as additional columns;\n\
  - run 'multiple_test_correction.py' on the p-values column afterward;\n\
  - basically, see here:\n\
https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.chi2_contingency.html\n"
  
synopsis2 = "detailed description:\n\
 0. Pre-requisite:\n\
  - numpy and scipy;\n\
 1. Input file and parameter:\n\
  - <input>: tab-delimited, with columns N to N+3 including integer values for;\n\
     2x2 contingency tables, one table per line; must have a header line;\n\
  - '-N N': the column where the 2x2 contingency table starts [1];\n\
  - will ;\n\
 2. Option(s), output, etc:\n\
  - '-g'|'--gtest': invoke 'lambda_=""log-likelihood""' option when running\n\
     'chi2_contingency()' from scipy [False];\n\
  - '-e'|'--scientific': p-values in scientific notations (e.g. 1.00e-2) [False];\n\
  - print results to stdout (lazy ...)\n\
  - will add options to accept larger contingency tables later (lazy ...)\n\
by ohdongha@gmail.com 20201110 ver 0.0\n\n"
#version_history
#20201110 ver 0.0 

#############################
### 0.2 parsing arguments ###
#############################
parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)
# positional arguments
parser.add_argument('input', type=argparse.FileType('r'))
# options
parser.add_argument('-N', dest="N", type=int, default= 1) 
parser.add_argument('-g', '--gtest', action="store_true", default=False) 
parser.add_argument('-e', '--scientific', action="store_true", default=False) 

args = parser.parse_args()

# parsing arguments
try:
	col_start = int( args.N ) - 1
	if col_start < 0:
		raise ValueError
except ValueError:
	print("'-N' expect a positive integer. Given an invalid input, default N=1 will be used.")
	col_start = 0
	
###########################
### 0.3 importing stuff ###
###########################
try:
	import numpy as np
	from scipy.stats import chi2_contingency
except ImportError as ErrorMessage:
	print(str(ErrorMessage))
	sys.exit()


###############
### 1. main ###
###############
num_line = 1
valid_line = True
header = True
for line in args.input:
	if header:
		print( "%s\tchi^2(g)\tp-value\tdof" % line.strip() )
		header = False
	else:
		tok = line.split('\t')
		num_line += 1
		try:
			obs = np.array( [ [ int(tok[col_start]), int(tok[col_start+1]) ], \
							  [ int(tok[col_start+2]), int(tok[col_start+3]) ] ] )
			if args.gtest:
				g, p, dof, expctd = chi2_contingency(obs, lambda_="log-likelihood")
			else:
				g, p, dof, expctd = chi2_contingency(obs)
			valid_line = True
		except ValueError:
			print("Line %d: non-integer values in columns %d to %d, skipping" % (num_line, col_start + 1, col_start + 4) )
			valid_line = False
			
		if valid_line:
			if args.scientific:
				print( "%s\t%f\t%.2e\t%d" % ( line.strip(), g, p, dof) )
			else:
				print( "%s\t%f\t%f\t%d" % ( line.strip(), g, p, dof) )
		else:
			print( "%s\tNA\tNA\tNA" % line.strip() )
