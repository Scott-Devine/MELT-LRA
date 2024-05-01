#!/usr/bin/env python3

# Compare combined and filtered freeze[4] VCF with results for individual samples
# and generate a freeze4-compatible MEI VCF.

import argparse
import csv
import gzip
import os
import re
import sys

# ------------------------------------------------------
# globals
# ------------------------------------------------------
VCF_REGEX = r'^.*\.vcf\.gz$'

# increase CSV max field size
csv.field_size_limit(256 * 1024 * 1024)

DEBUG = True

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
# read_individual_vcf
# ------------------------------------------------------
def read_individual_vcf(fpath, seq_index, ref_pos_index):
    lnum = 0
    debug("reading " + fpath)
    
    with gzip.open(fpath, 'rt') as fh:
        cr = csv.reader(fh, delimiter='\t')
        for row in cr:
            lnum += 1
            if re.match(r'^\#.*$', row[0]):
                continue

            # index each variant by ALT sequence and ref position
            (chrom, pos, vcf_id, ref, alt, qual, filt, inf, fmt, *rest) = row
            ins_seq = alt[1:]
            
            if ins_seq not in seq_index:
                seq_index[ins_seq] = []
            seq_index[ins_seq].append(row)

            ins_len = len(ins_seq)
            ref_pos = ":".join([chrom, pos, str(ins_len)])
            if ref_pos not in ref_pos_index:
                ref_pos_index[ref_pos] = []
            ref_pos_index[ref_pos].append(row)
            
# ------------------------------------------------------
# read_individual_vcf_dir
# ------------------------------------------------------
def read_individual_vcfs_dir(dpath):
    seq_index = {}
    ref_pos_index = {}
    nfiles = 0

    for file in os.listdir(dpath):
        if re.match(VCF_REGEX, file):
            fpath = os.path.join(dpath, file)
            read_individual_vcf(fpath, seq_index, ref_pos_index)
            nfiles += 1
#        if DEBUG and nfiles >= 20:
#            debug("reading only first 20 files")
#            break

    info("read " + str(nfiles) + " VCF file(s) from " + dpath)
    return { 'seq_index': seq_index, 'ref_pos_index': ref_pos_index }
    
# ------------------------------------------------------
# main()
# ------------------------------------------------------
def main():

    # input
    parser = argparse.ArgumentParser(description='Compare combined and filtered freeze[4] VCF with results for individual samples.')
    parser.add_argument('--individual_vcf_dir', required=True, help='Path to VCF files for individual samples/haplotypes.')
    parser.add_argument('--freeze_vcf', required=True, help='Path to merged and filtered VCF.')
    parser.add_argument('--output_vcf', required=True, help='Path to output MEI VCF based on merging freeze with individual samples.')
    parser.add_argument('--max_dist_bp', required=True, default=100, type=int, help='Maximum distance to merge an MEI with an exact sequence match to the reference.')
    args = parser.parse_args()

    # index individual VCF files
    vcf_ind = read_individual_vcfs_dir(args.individual_vcf_dir)
    rpi = vcf_ind['ref_pos_index']
    si = vcf_ind['seq_index']
    
    # iterate over merged freeze VCF
    info("reading " + args.freeze_vcf)

    # number of position matches with/without exact seq match
    exact_match = 0
    no_exact_match = 0

    # total number of expected matches based on genotypes ('1|1|.' = 2 matches expected)
    total_n_expected = 0
    # actual number of matches found
    total_n_found_ref_pos = 0
    total_n_found_ref_pos_seq = 0
    total_n_found_seq = 0
    
    # TODO - check that the matches are found in the haplotypes where we expect
    
    lnum = 0
    with gzip.open(args.freeze_vcf, 'rt') as fh:
        cr = csv.reader(fh, delimiter='\t')
        for row in cr:
            lnum += 1
            if re.match(r'^\#.*$', row[0]):
                continue

            (chrom, pos, vcf_id, ref, alt, qual, filt, inf, fmt, *rest) = row

            # examine insertions only
            if not re.match(r'.*-INS-.*', vcf_id):
                continue

            # TODO - expected haplotype count should be based on the number of MEIs that carry over
            # count '1's to find expected number of matches (not necessarily exact)
            n_expected = 0
            for gt in rest:
                for hap in gt.split('|'):
                    if hap == '1':
                        n_expected += 1
            total_n_expected += n_expected

            # -----------------------------------------------------------
            # 1. find matching MEIs (if any) by sequence position
            # -----------------------------------------------------------
            ins_len = len(alt) - 1
            ref_pos = ":".join([chrom, pos, str(ins_len)])
            exact_matches = {}
            
            if ref_pos in rpi:
                meis = rpi[ref_pos]
                nm = len(meis)
                total_n_found_ref_pos += nm
                
                # count how many have exact sequence match
                num_em = 0
                for mei in meis:
                    (m_chrom, m_pos, m_vcf_id, m_ref, m_alt, m_qual, m_flt, m_inf, m_fmt, *m_rest) = mei
                    if m_alt == alt:
                        num_em += 1
                        exact_matches[m_vcf_id] = True

                if num_em == 0:
                    no_exact_match += 1
                    warn("found " + str(nm) + " MEIs for " + ref_pos + " " + str(num_em) + " / " + str(nm) + " have exact ref position and ALT seq match")
                    for mei in meis:
                        warn(" mei=" + str(mei))
                    
                elif num_em > 0:
                    exact_match += 1

                total_n_found_ref_pos_seq += num_em

            # -----------------------------------------------------------
            # 2. find matching MEIs by sequence identity
            # -----------------------------------------------------------
            ins_seq = alt[1:]
            
            if ins_seq in si:
                meis = si[ins_seq]
                nm = len(meis)
#                info("found " + str(nm) + " MEIs with seq=" + alt)

                # filter out the exact matches that have already been found
                new_meis = []
                for mei in meis:
                    (m_chrom, m_pos, m_vcf_id, m_ref, m_alt, m_qual, m_flt, m_inf, m_fmt, *m_rest) = mei
                    if m_vcf_id not in exact_matches:
                        new_meis.append(mei)

                n_new = len(new_meis)
#                info("found " + str(n_new) + "/" + str(nm) + " new MEIs with seq=" + alt)

                # TODO- check how many of these are near enough to the original
                close_meis = []
                if n_new > 0:
#                    info("checking exact seq matches for proximity to " + vcf_id + " " + ref_pos)
                    for new_mei in new_meis:
                        (m_chrom, m_pos, m_vcf_id, m_ref, m_alt, m_qual, m_flt, m_inf, m_fmt, *m_rest) = new_mei
                        m_ins_len = len(m_alt) - 1
                        if chrom != m_chrom:
                            continue
                        dist = abs(int(pos) - int(m_pos))
                        if dist < args.max_dist_bp:
                            close_meis.append(new_mei)
                        
                            info(" found " + m_vcf_id + " " + m_chrom + " " + m_pos + " " + str(m_ins_len) + " dist=" + str(dist))
                        # TODO
                        
                n_close = len(close_meis)

                total_n_found_seq += n_close

            # TODO - what's left?
            #  -matches that are close in sequence and position but not identical in both
                
        info("exact seq match = " + str(exact_match))
        info("no exact seq match = " + str(no_exact_match)) 

        info("expected matches = " + str(total_n_expected))
        info("matches by reference position and sequence = " + str(total_n_found_ref_pos_seq))
        info("matches by reference position and insertion length = " + str(total_n_found_ref_pos))
        info("matches by sequence and proximity = " + str(total_n_found_seq))
                
        # TODO - select one to write to the output file
                
            # TODO - look for matches by sequence but NOT position and see how many are nearby

            # no matches found by reference position - check sequence index instead
#            else:
#                pass
            
#            if ref_pos not in ref_pos_index:
#                ref_pos_index[ref_pos] = []
#            ref_pos_index[ref_pos].append(row)
    

    # TODO - write output file
    
    # TODO
    # read merged/filtered VCF and compare:
    # for each position
    #   look up sample MEIs by position alone
    #   count how many MEIs match the sequence exactly
    #   print one (which one) of the MEIs with an exact sequence match

    #  count sample MEIs that match by sequence but not position and see how far away they are

    #  compare counts with number of 1s in the genotypes
    
    # for each position report/record:
    #   -actual number of exact occurrences (i.e., haplotype+sample count, with or without 'un')
    #   -number of merged occurrences (i.e., count up the '1's in the reported genotypes)
    
if __name__ == '__main__':
    main()
