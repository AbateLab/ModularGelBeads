#!/usr/bin/python
# encoding:utf8
# authors: Erik Garrison, Sébastien Boisvert
"""This script takes two fastq or fastq.gz files and interleaves them
Usage:
    interleave-fasta fasta_file1 fasta_file2
"""

import sys

def interleave(f1, f2, f_out):
    """Interleaves two (open) fastq files.
    """
    while True:
        line = f1.readline()
        if line.strip() == "":
            break
        f_out.write(line)
        
        for i in xrange(3):
             f_out.write(f1.readline())
        
        for i in xrange(4):
             f_out.write(f2.readline())

if __name__ == '__main__':
    try:
        file1 = sys.argv[1]
        file2 = sys.argv[2]
        file3 = sys.argv[3]
    except:
        print __doc__
        sys.exit(1)

    if file1[-2:] == "gz":
        import gzip
        with gzip.open(file1) as f1:
            with gzip.open(file2) as f2:
                with gzip.open(file3, 'wb') as f_out:
                    interleave(f1, f2, f_out)
    else:
        with open(file1) as f1:
            with open(file2) as f2:
                with open(file3, 'w') as f_out:
                    interleave(f1, f2, f_out)
       
