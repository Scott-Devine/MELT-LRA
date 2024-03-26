#!/usr/bin/env python3

# Check that the output from water alignments is complete.

import argparse
import csv
import os
import re
import sys

# ------------------------------------------------------
# globals
# ------------------------------------------------------
WATER_REGEX = r'^.*-water(-rev)?\.out$'
SAMPLE_REGEX = r'^((HG|NA)\d+).*-water(-rev)?\.out$'
FINAL_LINE_REGEX = r'^#---------------------------------------$'
DEBUG = False

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
# main()
# ------------------------------------------------------
def main():
    # input
    parser = argparse.ArgumentParser(description='Check that the output from water alignments is complete.')
    parser.add_argument('--water_dir', required=True, help='Path to directory containing water output.')
    args = parser.parse_args()

    if not os.path.exists(args.water_dir):
        fatal("specfiied --water_dir does not exist")
    if not os.path.isdir(args.water_dir):
        fatal("specfiied --water_dir is not a directory")

    # read water files in directory
    wfiles = [f for f in os.listdir(args.water_dir) if re.match(WATER_REGEX, f)]
    debug("files = " + str(wfiles))

    sample_files = {}

    # check that each sample has the same number of water output files
    for wfile in wfiles:
        wpath = os.path.join(args.water_dir, wfile)
        if not os.path.isfile(wpath):
            warn(wpath + " is not a file")
        m = re.match(SAMPLE_REGEX, wfile)
        if not m:
            fatal("couldn't parse sample id from water file " + wfile)
        sample = m.group(1)

        if sample not in sample_files:
            sample_files[sample] = []
        sample_files[sample].append(wpath)

    info("num_samples: " + str(len(sample_files)))
    n_files_per_sample = None
    
    for sample in sample_files:
        n_files = len(sample_files[sample])
        info(sample + ": " + str(n_files) + " files")
        if n_files_per_sample is None:
            n_files_per_sample = n_files
        elif files_per_sample != n_files:
            fatal("sample " + sample + " has " + str(n_files) + " output files, exepcted " + str(n_files_per_sample))
        
    # check that each contains the same number of alignments and ends cleanly
    n_alignments_per_file = None

    for wfile in wfiles:
        wpath = os.path.join(args.water_dir, wfile)

        # count alignments
        debug("checking " + wpath)
        n_alignments = 0
        
        # previous 2 lines
        prev_prev_line = None
        prev_line = None
        
        with open(wpath, 'r') as wfh:
            for line in wfh:
                if re.match(r'^# Aligned_sequences: 2$', line):
                    n_alignments = n_alignments + 1

                prev_prev_line = prev_line
                prev_line = line

        # check that final 2 lines match
        if not re.match(FINAL_LINE_REGEX, prev_prev_line) or not re.match(FINAL_LINE_REGEX, prev_line):
            fatal("file " + wpath + " appears to end prematurely")
                
        if n_alignments_per_file is None:
            n_alignments_per_file = n_alignments
        elif n_alignments != n_alignments_per_file:
            fatal("file " + wpath + " contains " + str(n_alignments) + " expected " + str(n_alignments_per_file))

    info("num_alignments_per_file: " + str(n_alignments_per_file))
            
if __name__ == '__main__':
    main()

