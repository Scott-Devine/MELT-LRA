#!/usr/bin/env python3

# Extract insertion sequences from PAV VCF and write them to FASTA.

import argparse
import csv
import gzip
import os
import re    
import sys

# increase CSV max field size
csv.field_size_limit(256 * 1024 * 1024)

def log_info(msg):
    sys.stderr.write("INFO - " + msg + "\n")
    sys.stderr.flush()

def log_fatal(msg):
    sys.stderr.write("FATAL - " + msg + "\n")
    sys.stderr.flush()
    sys.exit(1)

# ------------------------------------------------------
# extract_seqs
# ------------------------------------------------------

def extract_seqs(vcf, min_seqlen):
    n_seqs = 0
    n_extracted = 0
    
    with gzip.open(vcf, 'rt') as fh:
        cr = csv.reader(fh, delimiter='\t')
        for row in cr:
            if re.match(r'^#', row[0]):
                continue
            
            (chrom, pos, vcf_id, ref, alt, qual, filt, info, fmt, *rest) = row

            inf_d = {}
            for inf in info.split(';'):
                (k, v) = inf.split('=')
                if k in inf_d:
                    log_fatal("key " + k + " already seen in " + info)
                inf_d[k] = v

            if inf_d['SVTYPE'] != 'INS':
                continue

            # strip reference base from alt sequence
            n_seqs += 1
            
            if len(ref) != 1:
                log_fatal("len(REF) != 1")
            if alt[0] != ref:
                log_fatal("alt[0] (" + alt[0] + " != ref (" + ref + ")")

            alt = alt[1:]
            if len(alt) >= min_seqlen:
                n_extracted += 1
                print(">", vcf_id + " " + chrom + ":" + pos + " " + info)
                print(alt)
            
        log_info("extracted " + str(n_extracted) + "/" + str(n_seqs) + " sequences from insertion variants")
                
# ------------------------------------------------------
# main()
# ------------------------------------------------------

def main():

    # input
    parser = argparse.ArgumentParser(description='Extract insertion sequences.')
    parser.add_argument('--vcf', required=True, help='Path to input VCF file.')
    parser.add_argument('--min_seqlen', required=False, type=int, default=1, help='Minimum insertion sequence length.')
    args = parser.parse_args()
    extract_seqs(args.vcf, args.min_seqlen)
        
if __name__ == '__main__':
    main()



