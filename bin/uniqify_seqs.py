#!/usr/bin/env python3

# Remove duplicate sequences from a set of FASTA files.

import argparse
import csv
import gzip
import os
import re
import sys

# ------------------------------------------------------
# globals
# ------------------------------------------------------
FASTA_SUFFIX_RE = r'\.(fa.gz|fsa.gz|fasta.gz)$'
FASTA_FILE_RE = r'^(.*)' + FASTA_SUFFIX_RE

# ------------------------------------------------------
# logging
# ------------------------------------------------------
def fatal(msg):
    sys.stderr.write("FATAL - " + msg + "\n")
    sys.exit(1)

def info(msg):
    sys.stderr.write("INFO - " + msg + "\n")

def debug(msg):
    if DEBUG:
        sys.stderr.write("DEBUG - " + msg + "\n")

def warn(msg):
    sys.stderr.write("WARN - " + msg + "\n")
    
# ------------------------------------------------------
# read_fasta_file
# ------------------------------------------------------
def read_fasta_file(fpath):
    seqs = []
    defline = None
    seq = ''

    def process_seq():
        nonlocal defline, seq
        seqs.append({ 'defline': defline, 'seq': seq, 'len': len(seq) })
        defline = None
        seq= ''
    
    with gzip.open(fpath, 'rt') as fh:
        for line in fh:
            # defline
            m = re.match(r'^>(.*)$', line)
            if m:
                if defline is not None:
                    process_seq()
                defline = m.group(1)
            else:
                rsl = line.rstrip()
                seq = seq + rsl
                rsl = rsl.upper()
                rsl = re.sub(r'[\s\-]+', '', rsl)

    if defline is not None:
        process_seq()

    info("read " + str(len(seqs)) + " sequences from " + fpath)
    return seqs

# ------------------------------------------------------
# read_fasta_dir
# ------------------------------------------------------
def read_fasta_dir(dpath):
    seqs = []

    for file in os.listdir(dpath):
        fpath = os.path.join(dpath, file)
        fseqs = read_fasta_file(fpath)
        seqs.extend(fseqs)

    return seqs
    
# ------------------------------------------------------
# main()
# ------------------------------------------------------
def main():
    # input
    parser = argparse.ArgumentParser(description='Remove duplicate sequences from a set of FASTA files.')
    parser.add_argument('--file_list', required=True, help='Path to flat file containing list of FASTA files or FASTA file directories to include.')
    args = parser.parse_args()

    # dictionary of all sequences
    seq_d = {}
    n_seqs = 0
    n_unique_seqs = 0

    def process_seq(s):
        nonlocal seq_d, n_seqs
        if s['seq'] in seq_d:
            seq_d[s['seq']].append(s['defline'])
            pass
        else:
            seq_d[s['seq']] = [s['defline']]
        n_seqs += 1
        
    # read list of files and/or directories
    with open(args.file_list) as flfh:
        for line in flfh:
            path = line.strip()
            seqs = []
            
            if not os.path.exists(path):
                fatal("file/directory does not exist: " + path)
            if os.path.isdir(path):
                seqs = read_fasta_dir(path, )
            else:
                seqs = read_fasta_file(path)

            for seq in seqs:
                process_seq(seq)
                
    info("n_seqs=" + str(n_seqs))
    info("n_unique_seqs=" + str(len(seq_d)))
    
if __name__ == '__main__':
    main()
