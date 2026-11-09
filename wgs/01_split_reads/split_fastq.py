#!/usr/bin/env python3
import gzip
from subprocess import Popen
import sys
import os

stem = sys.argv[1]
directory = sys.argv[2]
fileout_stem = sys.argv[3]
nb_reads = int(sys.argv[4])

file1_path_r1 = os.path.join(directory, stem + '.R1.fastq.gz')
file2_path_r2 = os.path.join(directory, stem + '.R2.fastq.gz')
file1_path_r1_alt = os.path.join(directory, stem + '_1.fq.gz')
file2_path_r2_alt = os.path.join(directory, stem + '_2.fq.gz')

try:
    file1 = gzip.open(file1_path_r1, 'rt')
    file2 = gzip.open(file2_path_r2, 'rt')
except FileNotFoundError:
    file1 = gzip.open(file1_path_r1_alt, 'rt')
    file2 = gzip.open(file2_path_r2_alt, 'rt')

count_read = 0
count_file = 1

head1 = file1.readline()
head2 = file2.readline()

fileout1 = open(f'{fileout_stem}_split_{count_file}_R1.fq', 'w')
fileout2 = open(f'{fileout_stem}_split_{count_file}_R2.fq', 'w')

while head1:
    seq1 = file1.readline()
    seq2 = file2.readline()
    line1 = file1.readline()
    line2 = file2.readline()
    qual1 = file1.readline()
    qual2 = file2.readline()

    fileout1.write(head1 + seq1 + line1 + qual1)
    fileout2.write(head2 + seq2 + line2 + qual2)
    count_read += 1

    head1 = file1.readline()
    head2 = file2.readline()

    if count_read == nb_reads:
        fileout1.close()
        fileout2.close()
        Popen(f'gzip {fileout_stem}_split_{count_file}_R1.fq', shell=True).wait()
        Popen(f'gzip {fileout_stem}_split_{count_file}_R2.fq', shell=True).wait()

        count_file += 1
        fileout1 = open(f'{fileout_stem}_split_{count_file}_R1.fq', 'w')
        fileout2 = open(f'{fileout_stem}_split_{count_file}_R2.fq', 'w')
        count_read = 0

fileout1.close()
fileout2.close()
Popen(f'gzip {fileout_stem}_split_{count_file}_R1.fq', shell=True).wait()
Popen(f'gzip {fileout_stem}_split_{count_file}_R2.fq', shell=True).wait()
file1.close()
file2.close()
