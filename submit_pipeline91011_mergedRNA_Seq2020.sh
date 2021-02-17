#!/bin/bash
#SBATCH -J Pipeline
#SBATCH -c 24
#SBATCH --mem=100000

# <script> <$2 inputfile> <$3 processedfile> <$3 resultsfile>


#####################
# Running scripts
#####################
# Arabidopsis
/Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/scripts/RNA_Seq2020.sh /Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/rawfiles/RNA_Seq/At_mergedraw_files /Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/processed/RNA_At /Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/results/results_At

# E. salsugineum- JGI
/Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/scripts/RNA_Seq2020.sh /Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/rawfiles/RNA_Seq/EsJ_mergedraw_files /Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/processed/RNA_Es /Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/results/results_Es


# Sisymbrium irio
/Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/scripts/RNA_Seq2020.sh /Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/rawfiles/RNA_Seq/Si_mergedraw_files /Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/processed/RNA_Si /Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/results/results_Si

# S. parvula
/Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/scripts/RNA_Seq2020.sh /Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/rawfiles/RNA_Seq/Sp_mergedraw_files /Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/processed/RNA_Sp /Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/results/results_Sp

