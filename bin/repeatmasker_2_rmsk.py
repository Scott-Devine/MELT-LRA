#!/usr/bin/env python3

# Convert raw RepeatMasker output to UCSC table browser format (approximately.)

import argparse
import csv
import gzip
import os
import re    
import sys

# increase CSV max field size
csv.field_size_limit(256 * 1024 * 1024)

def info(msg):
    sys.stderr.write("INFO - " + msg + "\n")
    sys.stderr.flush()

def fatal(msg):
    sys.stderr.write("FATAL - " + msg + "\n")
    sys.stderr.flush()
    sys.exit(1)

    
# ------------------------------------------------------
# main()
# ------------------------------------------------------

def main():

    # input
    parser = argparse.ArgumentParser(description='Convert raw RepeatMasker output to UCSC table browser format.')
    parser.add_argument('--repeatmasker', required=True, help='Path to to gzipped RepeatMasker output file.')
    args = parser.parse_args()

    fr = re.compile(r'\s+')
    lnum = 0
    
    with gzip.open(args.repeatmasker, 'rt') as rfh:
        for line in rfh:
            lnum = lnum + 1
            if re.match(r'^\s*((SW |score ).*)?$', line):
                continue
            line = line.strip()
            fields = fr.split(line)
            nf = len(fields)
            if nf != 15:
                fatal("wrong number of fields at line " + str(lnum) + ": " + str(nf))

            # strand is either '+' or 'C'
            [sw_score, pct_div, pct_del, pct_ins, seq, begin, end, left, strand, repeat, rclass, rep_begin, rep_end, rep_left, rep_id] = fields

            # parse repeat class/family
            rep_class = rclass
            rep_family = rclass
            m = re.match(r'^(.*)\/(.*)$', rclass)
            if m:
                rep_class = m.group(1)
                rep_family = m.group(2)

            # output
            # UCSC bin / not used
            ucsc_bin = 1
            milliDiv = float(pct_div) * 10.0
            milliDel = float(pct_del) * 10.0
            milliIns = float(pct_ins) * 10.0
            genoName = seq
            genoStart = int(begin) - 1
            genoEnd = end
            genoLeft = left
            ucsc_strand = '+' if strand == '+' else '-'
            
            ocols = [ucsc_bin, sw_score, milliDiv, milliDel, milliIns, genoName, genoStart, genoEnd, genoLeft, ucsc_strand, repeat, rep_class, rep_family, rep_begin, rep_end, rep_left, rep_id]
            print("\t".join([str(oc) for oc in ocols]))

    info("read " + str(lnum) + " line(s) from " + args.repeatmasker)


if __name__ == '__main__':
    main()



