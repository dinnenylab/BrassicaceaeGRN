#!/usr/bin/env python
import sys, os, argparse
from argparse import RawTextHelpFormatter

##############################
### 0.1 script description ###
##############################
synopsis1 = "\
  - the main script for PiP (Phylogenetically informed Profiling);\n\
  - draw pairs plot (seaborn.pairplot) for a matrix of numbers;\n\
  - print Peason and Spearman correlation matrices for all pairs;\n\
  - can iterate for multiple subsets of the input matrix;\n"
  
synopsis2 = "detailed description:\n\
 0. Pre-requisite:\n\
  - python3, seaborn 0.10, panda, and matplotlib;\n\
 1. Input files:\n\
  - <input>: tab-delimited, with the 1st column including geneID (or equivalent,\n\
     e.g. ortholog pair ID), followed by numerical columns; a header with column\n\
     names (e.g. sampeID) required;\n\
  - with '-N N', assumes columns 1 to N-1 containing geneID (or equivalent) to\n\
     identify rows to be included in subsets (given by '-s'); useful when using\n\
     pairwise orthologs for calculating correlations among different species;\n\
  - <subset> (optional, given by '-s') is tab-delimited, with the 1st column\n\
     subset names (e.g. GO terms), followed by row names (e.g. geneID) in each\n\
     subset, separated by '|' (i.e. as in a BiNGO output file); can have multiple\n\
     lines with the same subset names.\n\
  *** note: when subsetting, all geneIDs are converted to UPPERCASE ... this is\n\
     because .bgo somehow print gene names in uppercase when asked to print all...\n\
 2. Options and parameters:\n\
  - '-s subset': see 'Input files:'; if not specified, work only once with all\n\
     rows in <input>;\n\
  - '-N N': the data to calculate correlations starts from Nth column, and\n\
     all values in 1st to (N-1)th columns are considered to select subsets [2];\n\
  - '-m min_rows_subset': print 'NA' for subsets with members less than [10];\n\
  - '-o output_name': if not specified, <output_name> = <input> w/o extension\n\
     + '_pplt.png';\n\
  - '-1'|'--one_liner': '-1' is the number 'one' NOT the lowercase of 'L';\n\
     correlation matrix as one-liner, also print p-values using scipy [False];\n\
  - '-n'|'--no_graph': skip pairs plot and print only correlation matrix [False];\n\
  - '-p'|'--print_subset': print a matrix of presence (=1) and absence (=0) of\n\
     rows for each subset as <input> w/o extension + '_tags.txt'; does not\n\
     calcurate correlations nor draw pairplots [False];\n\
  - '-c'|'--draw_CI': add confidence intervals (95%) to the pair plot [False];\n\
  - '-xylim': values to use for xlim=ylim; use '-xylim=-1' to draw pairplots\n\
     with flexible x, y limits (useful for detecting outliers) [2.25];\n\
  - '--SVG': print pairplots in .sgv format, instead of .png [False];\n\
# - '-pplt_opt': string of options to be passed to seaborn.pairplot.set();\n\
#    put in quotation marks; e.g. -pplt_opt 'xlim=(-2, 2), ylim=(-2, 2)'\n\
 4. Output:\n\
  - print the pairs plot as <output_name>_pplt.png\n\
  - print Pearson and Spearman correlation matrices to stdout\n\
by ohdongha@gmail.com 2020929 ver 0.2.3\n\n"
#version_history
#20201114 code renamed to "PiP_correlation_matrix_pairwise_orthologs.py" PiP for "Phylogenetically informed Profiling"
#20200929 ver 0.2.3 pplt_code has a line to control font sizes
#20200905 ver 0.2.2 fixed '-p'|'--print_subset' option to deal with pairwise orthologs input + SVG output
#20200831 ver 0.2.1 added options not to set xlim and ylim (useful when detecting outliers) + to add CI trendline, etc
#20200823 ver 0.2 modified to deal with pairwise orthologs input (rather than an all 1-vs-1 matrix)
#20200611 ver 0.1.4 accept xlim = ylim value as argument ...
#20200608 ver 0.1.3 print also p-values for Pearson and Spearman correlations when using '-1' option ...
#20200602 ver 0.1.2 use sns.PairGrid instead of sns.pairplot to enable -xlim, -ylim, etc ...
#20200530 ver 0.1.1 '-p'|'--print_subset' option added
#20200529 ver 0.1 '-s subset' option added
#20200528 ver 0.0.2 '-1'|'--one_liner' option added
#20200527 ver 0.0.1 now using python 3.8 and seaborn 0.10+
#20200526 ver 0.0 modified from plot_hist_dist_2D.py 

#############################
### 0.2 parsing arguments ###
#############################
parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)
# positional arguments
parser.add_argument('input', type=argparse.FileType('r'))
#parser.add_argument('col1', type=int)
#parser.add_argument('col2', type=int)
# options
parser.add_argument('-s', dest="subset", type=str, default= "__NA__")
parser.add_argument('-N', dest="N", type=int, default= 2) 
parser.add_argument('-m', dest="min_rows_subset", type=int, default= 10) 
parser.add_argument('-o', dest="output_name", type=str, default= "__NA__") 
parser.add_argument('-1', '--one_liner', action="store_true", default=False) 
parser.add_argument('-n', '--no_graph', action="store_true", default=False) 
parser.add_argument('-p', '--print_subset', action="store_true", default=False) 
parser.add_argument('-c', '--draw_CI', action="store_true", default=False) 
parser.add_argument('--SVG', action="store_true", default=False) 
parser.add_argument('-xylim', dest="xylim", type=float, default= 2.25) 
#parser.add_argument('-pplt_opt', dest="pplt_opt", type=str, default= "__NA__") 

args = parser.parse_args()

# parsing arguments
if args.output_name == "__NA__":
	outfile_base = os.path.splitext(args.input.name)[0] # add "_pplt.png" etc later
else:
	outfile_base = args.output_name

if args.subset == "__NA__":
	work_with_subset = False
else:
	work_with_subset = True

try:
	xylim = float( args.xylim )
except ValueError:
	print("invalid value for -xylim, using default=2.25")
	xylim = 2.25

try:
	col_start = int( args.N ) - 1
	if col_start < 1:
		raise ValueError
except ValueError:
	print("'-N' expect an integer larger than 1. Given an invalid input, default N=2 will be used.")
	col_start = 1
	
	
###########################
### 0.3 importing stuff ###
###########################
try:
	import numpy as np, pandas as pd
	import seaborn as sns; sns.set(style="ticks", color_codes=True)
	import matplotlib.pyplot as plt
	plt.switch_backend('agg') 
	import scipy.stats as ss

except ImportError as ErrorMessage:
	print(str(ErrorMessage))
	sys.exit()

def hide_current_axis(*args, **kwds):
	plt.gca().set_visible(False)


##############################################
### 1. read input (and subset) and iterate ###
##############################################
## 1.1 reading input 
df=pd.read_csv(args.input,delimiter='\t',header=0, low_memory=False) # let's think about dtype issues later ...
#print(list(df.columns.values)) #file header
#df.describe()

## 1.2 reading subset
if work_with_subset:
	try:
		fin_subset = open( args.subset, "r")
	except FileNotFoundError:
		print( "## Subset file: %s not found, working with the entire set" % args.subset )
		work_with_subset = False
		
if work_with_subset:
	subset_dict = dict() # key = subset_name; value = set of members (geneIDs)
	subset_list = [] # list of subset names as they appear in the subset file
	for line in fin_subset:
		tok = line.split('\t')
		if tok[0] not in subset_dict:
			subset_dict[ tok[0] ] = set( tok[1].upper().split('|') ) # convert all to uppercase -_-;;; 
			subset_list.append( tok[0] )
		else:
			subset_dict[ tok[0] ] = subset_dict[ tok[0] ] | set( tok[1].upper().split('|') )
	print( "## Finished reading %d subsets from %s"  % (len( subset_dict ), fin_subset.name ) )
	fin_subset.close()


##################################################
### 2. print pairs plot and correlation matrix ###
##################################################
#current sns.PairGrid script to print pairs plot - modify as needed
#v0.2 added "dropna=True" ... will this work as intended? ... yes it seems like
#v0.2.3 added "g = sns.set(font_scale=2);\" to increase tick font sizes
if xylim < 0:
	pplt_code = "g = sns.PairGrid(df_subset,dropna=True); \
					g = g.map_lower(plt.scatter,color='.3',s=14,alpha=0.2); \
					g = g.map_upper(hide_current_axis); \
					g = g.map_diag(sns.kdeplot,color='.3',shade=True)"
else: 
	pplt_code = "g = sns.set(font_scale=2); \
					g = sns.set_style('ticks'); \
					g = sns.PairGrid(df_subset,dropna=True); \
					g = g.map_lower(plt.scatter,color='.3',s=14,alpha=0.2); \
					g = g.map_upper(hide_current_axis); \
					g = g.map_diag(sns.kdeplot,color='.3',shade=True); \
					g = g.set(xlim=(-%f,%f),ylim=(-%f,%f))" % (xylim, xylim, xylim, xylim)
if args.draw_CI:
	pplt_code=pplt_code.replace("(plt.scatter,color='.3',s=14,alpha=0.2)", \
                                "(sns.regplot,ci=95,fit_reg=True,color='.3',scatter_kws={'s':14,'alpha':0.2})")

if args.SVG:
	savefig_code = "g.savefig(outfile, format='svg')"
else:
	savefig_code = "g.savefig(outfile, dpi = 300)"

#pplt_code = "g = sns.PairGrid(df_subset); \
#				g = g.map_upper(plt.scatter, s=17, alpha=0.3); \
#				g = g.map_lower(sns.kdeplot, cmap='magma'); \
#				g = g.map_diag(sns.kdeplot, shade=True); \
#				g = g.set(xlim=(-2.25,2.25), ylim=(-2.25,2.25))"
#layout_code = "plt.tight_layout()" # seems not working

if work_with_subset:
	start = True
	if args.print_subset:
		presence_in_subset_dict = dict() # k = subset_name; v = dict with k = geneID (@1stCol); v = 0 or 1 
		geneID_list = df[ df.columns[0] ].tolist() # list of geneIDs @1stCol, use these as anchors			
		
	for subset_name in subset_dict:
		for n in range(0, col_start): # v0.2 search multiple columns to subset; default col_start=1 (i.e. 2nd column, N=2)
			if n == 0:
				df_subset = df.loc[ df[df.columns[0]].str.upper().isin( subset_dict[subset_name] ) ]
			else:
				df_subset = pd.concat( [ df_subset, df.loc[ df[df.columns[n]].str.upper().isin( subset_dict[subset_name] ) ] ] ).drop_duplicates().reset_index(drop=True)
				
		num_rows = df_subset.shape[0]
		if num_rows >= args.min_rows_subset: # min_rows_subset is not very useful now ... whatever ...
			num_pairs = 0 # this will be more useful

			## with '--print_subset', remember presence/absence of each geneID in each subset (=GO term)
			## and skip everything else ... v0.2.2
			if args.print_subset:
				presence_in_subset_dict[subset_name] = dict() # k = geneID (@1stCol); v = 0 or 1
				num_members = 0
				for g in geneID_list:
					if g in df_subset[ df_subset.columns[0] ].tolist(): # v0.2.2
						presence_in_subset_dict[subset_name][g] = "1"
						num_members += 1
					else:
						presence_in_subset_dict[subset_name][g] = "0"
				print("found %d rows (=gene or ortholog group/pair IDs) in subset %s" % (num_members, subset_name) ) 
			else:
				## print pairs plot (unless told not to)
				if not args.no_graph:
					if args.SVG:
						outfile = "pplt_%s_" % subset_name + outfile_base + ".svg" 
					else:
						outfile = "pplt_%s_" % subset_name + outfile_base + ".png" 
					exec( pplt_code )
	#				exec( layout_code )		
					exec( savefig_code )	
					plt.close()
	
				## calculate and print correlation matrix
				if args.one_liner:
					## calculate and print correlation and p-values using scipy
					col_list = list(df_subset)[col_start:] # data start from Nth column
					num_col = len(col_list) 
					ij_list = [ (i,j) for i in range(0, num_col) for j in range(i+1, num_col) ] # to iterate
					
					df_subset_array_dict = dict() # key = colname, value = numpy.ndarray of the column values
					N_list = [] # list of number of non-NaN pairs for each comparison 
					Pr_list = [] # list of Pcorr coefficients 
					Pp_list = [] # list of Pcorr p-values 
					Sr_list = [] # list of Pcorr coefficients 
					Sp_list = [] # list of Pcorr p-values 
					
					for col in col_list:
						df_subset_array_dict[ col ] = df_subset[ col ].to_numpy()
					for i, j in ij_list:
						x = df_subset_array_dict[ col_list[i] ]
						y = df_subset_array_dict[ col_list[j] ]
						bad = ~np.logical_or(np.isnan(x), np.isnan(y)) # to deal with NaN (https://stackoverflow.com/a/48591908/6283377)
						num_pairs = np.count_nonzero(bad == True) # count pairs with no NaN
						N_list.append( num_pairs ) 
	
						if num_pairs >= args.min_rows_subset:
							Pr, Pp = ss.pearsonr( np.compress(bad, x), np.compress(bad, y) ) 
							Sr, Sp = ss.spearmanr( np.compress(bad, x), np.compress(bad, y) )
							Pr_list.append("%.4f" % Pr)
							Pp_list.append("%.2E" % Pp)
							Sr_list.append("%.4f" % Sr)
							Sp_list.append("%.2E" % Sp)
						else:
							Pr_list.append("NA")
							Pp_list.append("NA")
							Sr_list.append("NA")
							Sp_list.append("NA")						
						
					## print correlations and p-values as one-liners		
					if start: # print header only once at the beginning
						print( "data\tsubset\tvalue", end='' ) # print a header			
						for i, j in ij_list:
							print( "\t%s-%s" % ( col_list[i], col_list[j] ), end='' )
						print()
						start=False
					print( "%s\t%s\t#ValidPairs\t%s" % (outfile_base, subset_name, '\t'.join(str(n) for n in N_list) ) ) 
					print( "%s\t%s\tPcorr\t%s" % (outfile_base, subset_name, '\t'.join(Pr_list) ) ) 
					print( "%s\t%s\tPcorr.p\t%s" % (outfile_base, subset_name, '\t'.join(Pp_list) ) ) 
					print( "%s\t%s\tScorr\t%s" % (outfile_base, subset_name, '\t'.join(Sr_list) ) ) 
					print( "%s\t%s\tScorr.p\t%s" % (outfile_base, subset_name, '\t'.join(Sp_list) ) ) 

				## without '--one_liner', print correlation as matrices (no p-values) 			
				else:
					df_subset = df_subset[ df.columns[col_start:] ]
					PcorrMatrix = df_subset.corr(method='pearson')
					ScorrMatrix = df_subset.corr(method='spearman')
					print( "\ndata=%s, subset=%s, #rows=%d" % (outfile_base, subset_name, num_rows ) )
					print( "Pearson:" )
					print( PcorrMatrix )
					print( "Spearman:" )
					print( ScorrMatrix )
				
			
	## after iterating over all subsets,
	## print presence/absence of each row (gene) in each subset (GO term)
	if args.print_subset:
		fout_tags = open( outfile_base + ".tags.txt", 'w')	
		
		header_list = [df.columns[0].strip()]
		for k in subset_list:
			if k in presence_in_subset_dict:
				header_list.append( k )
		fout_tags.write( '\t'.join( header_list ) + '\n') # wrting the header
	
		for g in geneID_list:
			presence_list = [g]
			for k in subset_list:
				if k in presence_in_subset_dict:
					presence_list.append( presence_in_subset_dict[k][g] )
			fout_tags.write( '\t'.join( presence_list ) + '\n') # writing presence/absence
	
		sum_list = ["sum"]
		for k in subset_list:
			if k in presence_in_subset_dict:
				sum_list.append( str( sum(x == "1" for x in presence_in_subset_dict[k].values()) ) )
		fout_tags.write( '\t'.join( sum_list ) + '\n') # writing the sum of all rows (genes) with "1" for each subset (GO term)
	
		fout_tags.close()
	
else: # without subset
	df_subset = df[ df.columns[col_start:] ]
	num_rows = df_subset.shape[0]

	## print pairs plot (unless told not to)
	if not args.no_graph:
		outfile = outfile_base + "_pplt.png"
		exec( pplt_code )					
		exec( savefig_code )	
		plt.close()

	## calculate and print correlation matrix
	if args.one_liner:
		## calculate Pearson correlation and p-values using scipy
		col_list = list(df_subset)
		num_col = len(col_list) 
		ij_list = [ (i,j) for i in range(0, num_col) for j in range(i+1, num_col) ] # to iterate
		
		df_subset_array_dict = dict() # key = colname, value = numpy.ndarray of the column values
		N_list = [] # list of number of non-NaN pairs for each comparison 
		Pr_list = [] # list of Pcorr coefficients 
		Pp_list = [] # list of Pcorr p-values 
		Sr_list = [] # list of Pcorr coefficients 
		Sp_list = [] # list of Pcorr p-values 
		
		for col in col_list:
			df_subset_array_dict[ col ] = df_subset[ col ].to_numpy()
		for i, j in ij_list:
			x = df_subset_array_dict[ col_list[i] ]
			y = df_subset_array_dict[ col_list[j] ]
			bad = ~np.logical_or(np.isnan(x), np.isnan(y)) # to deal with NaN (https://stackoverflow.com/a/48591908/6283377)
			Pr, Pp = ss.pearsonr( np.compress(bad, x), np.compress(bad, y) ) 
			Sr, Sp = ss.spearmanr( np.compress(bad, x), np.compress(bad, y) )
#			Pr, Pp = ss.pearsonr( x, y ) 
#			Sr, Sp = ss.spearmanr( x, y )
			N_list.append( np.count_nonzero(bad == True) ) # count pairs with no NaN
			Pr_list.append("%.4f" % Pr)
			Pp_list.append("%.2E" % Pp)
			Sr_list.append("%.4f" % Sr)
			Sp_list.append("%.2E" % Sp)
			
		## print Pearson and Spearman correlations and p-values as one-liners		
		print( "data\tvalue", end='' ) # print a header			
		for i, j in ij_list:
			print( "\t%s-%s" % ( col_list[i], col_list[j] ), end='' )
		print( "\n%s\t#ValidPairs\t%s" % (outfile_base, '\t'.join(str(n) for n in N_list) ) ) 
		print( "%s\tPcorr\t%s" % (outfile_base, '\t'.join(Pr_list) ) ) 
		print( "%s\tPcorr.p\t%s" % (outfile_base, '\t'.join(Pp_list) ) ) 
		print( "%s\tScorr\t%s" % (outfile_base, '\t'.join(Sr_list) ) ) 
		print( "%s\tScorr.p\t%s" % (outfile_base, '\t'.join(Sp_list) ) ) 
				
	else: # if not one-liner, just print the easy way ... it seems automatically deals with NaN? 
		PcorrMatrix = df_subset.corr(method='pearson')
		ScorrMatrix = df_subset.corr(method='spearman')
		print( "\ninput=%s, #rows=%d" % (args.input.name, num_rows ) )
		print( "Pearson:" )
		print( PcorrMatrix )
		print( "Spearman:" )
		print( ScorrMatrix )
