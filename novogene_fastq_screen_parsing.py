#!/usr/bin/env python3

import os

data_dir = '/user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/results/fastq_screen/'

out_file = f"{data_dir}/fastq_screen_summary.csv"

out_fh = open(out_file, 'w')
out_fh.write('organism,sample,type,percentage\n')

for text_file in os.listdir(data_dir):
    if not text_file.endswith(".txt"):
        continue
    
    sample = text_file.replace("_screen.txt", "")
    text_file = f"{data_dir}/{text_file}"
    screen_fh = open(text_file, 'r')
    for line in screen_fh:
        line = line.rstrip()
        data = line.split("\t")
        if line.startswith("Mouse") or line.startswith("Tbru"):
            organism = data[0]
            print(organism)
            out_fh.write(f"{organism},{sample},01_unmapped_pc,{data[3]}\n")
            out_fh.write(f"{organism},{sample},05_one_hit_one_genome,{data[5]}\n")
            out_fh.write(f"{organism},{sample},04_multiple_hits_one_genome,{data[7]}\n")
            out_fh.write(f"{organism},{sample},03_one_hit_multiple_genomes,{data[9]}\n")
            out_fh.write(f"{organism},{sample},02_multiple_hits_multiple_genomes,{data[11]}\n")

out_fh.close()
