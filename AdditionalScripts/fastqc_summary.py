#!/usr/local/bin/python
import sys

# USAGE: python fastqc_summary.py sample_names.txt_dir outdir
sample_names_dir = sys.argv[1].replace('\\','')
outdir = sys.argv[2]

f = open(sample_names_dir, 'r')
sample_names = []
for file in f:
    sample_names.append(file.strip())

out = open(outdir, 'w')
for sample in sample_names:
    sample_sum = open( sample + '/summary.txt', 'r')
    for line in sample_sum:
        line = line.strip()
        out.write(line + '\n')

