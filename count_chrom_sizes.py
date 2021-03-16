#!/usr/bin/env python
#count_chrom_sizes.py input.fa output
import sys

fin = open(sys.argv[1], 'r')
#fout = open(sys.argv[2], 'w')
#fout2 = open(sys.argv[2]+".list", 'w')
fout2 = open(sys.argv[2], 'w')

seq_size = 0
head = ''
seq = list()
DNAbases = set('ATGCNatgcnRYSWKMBDHVryswkmbdhv')

for line in fin:
    if line[0] == '>':     # beginning of a sequence
        # before processing, the last sequence must be written
        if(head != ''):
#           new_head = '>' + head[:-1] + '\t' + str(len(seq)) + '\n' 
            new_head = head[:-1] + '\t' + str(len(seq)) + '\n' 
#           fout.write(new_head)
            fout2.write(new_head)
#           fout.write("".join(seq)+ '\n')
        # start a new sequence
        seq = list()
        head = line[1:]
    else: # in the middle of sequence
        for c in line:
            if(c in DNAbases):
                seq.append(c)

# write the last seq
#new_head = '>' + head[:-1] + '\t' + str(len(seq)) + '\n' 
new_head = head[:-1] + '\t' + str(len(seq)) + '\n' 
#fout.write(new_head)
fout2.write(new_head)
#fout.write("".join(seq)+ '\n')


