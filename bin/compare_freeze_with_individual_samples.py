#!/usr/bin/env python3

# Compare combined and filtered freeze[4] VCF with results for individual samples
# and generate a freeze4-compatible MEI VCF.

import argparse
import csv
from datetime import date
import gzip
import os
import re
import sys

# ------------------------------------------------------
# globals
# ------------------------------------------------------
VERSION = '1.4.2'

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

    m = re.match(r'.*((NA|HG)\d+)-(h1|h2|un).*', fpath)
    if not m:
        fatal("couldn't parse sample/haplotype from " + fpath)

    sample = m.group(1)
    haplotype = m.group(3)
    info_lines = []
    
    with gzip.open(fpath, 'rt') as fh:
        cr = csv.reader(fh, delimiter='\t')
        for row in cr:
            lnum += 1
            if re.match(r'^\#.*$', row[0]):
                if re.match(r'^\#\#INFO=.*$', row[0]):
                    info_lines.append(row[0])
                continue

            # index each variant by ALT sequence and ref position
            (chrom, pos, vcf_id, ref, alt, qual, filt, inf, fmt, *rest) = row

            # add sample/haplotype to vcf_id
            row[2] = ":".join([sample, haplotype, vcf_id])

            ins_seq = alt[1:]
            
            if ins_seq not in seq_index:
                seq_index[ins_seq] = []
            seq_index[ins_seq].append(row)

            ins_len = len(ins_seq)
            ref_pos = ":".join([chrom, pos, str(ins_len)])
            if ref_pos not in ref_pos_index:
                ref_pos_index[ref_pos] = []
            ref_pos_index[ref_pos].append(row)

    return { 'info': info_lines }
            
# ------------------------------------------------------
# read_individual_vcf_dir
# ------------------------------------------------------
def read_individual_vcfs_dir(dpath):
    seq_index = {}
    ref_pos_index = {}
    nfiles = 0
    info_lines = None
    
    for file in os.listdir(dpath):
        if re.match(VCF_REGEX, file):
            fpath = os.path.join(dpath, file)
            vc = read_individual_vcf(fpath, seq_index, ref_pos_index)
            if info_lines is None:
                info_lines = vc['info']
            nfiles += 1

    info("read " + str(nfiles) + " VCF file(s) from " + dpath)
    return { 'seq_index': seq_index, 'ref_pos_index': ref_pos_index, 'info': info_lines }
    
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
    info_lines = vcf_ind['info']
    
    # iterate over merged freeze VCF
    info("reading " + args.freeze_vcf)

    # output file
    ofh = open(args.output_vcf, 'w')
    info("writing " + args.output_vcf)
    ofh.write("##fileformat=VCFv4.2\n")
    ofh.write("##fileDate=" + date.today().strftime("%Y%m%d") + "\n")
    ofh.write("##source=MELT-RISC " + VERSION + "\n")
    ofh.write("##reference=" + args.freeze_vcf + "\n")
    
    # number of freeze VCF positions matched (including insertion length)
    n_freeze_vcf_pos_plus_inslen_found = 0
    # number of freeze VCF insertion seqs matched
    n_freeze_vcf_ins_seqs_found = 0
    # reference positions found, irrespective of insertion length or sequence
    freeze_vcf_pos_found = {}
    # ME types found
    ME_types_found = { 'ALU': 0, 'SVA': 0, 'LINE1': 0 }
    
    # counts by reference entry (ref pos including insertion length)
    entry_counts = {
        'exact': 0,
        'exact_seq_close_pos' : 0,
        'exact_pos_close_seq': 0,
        'total': 0,                 # total based on entries with either ref position or ref seq found in index
    }
    
    # counts by reference haplotype
    hap_counts = {
        'exact': 0,
        'exact_seq_close_pos' : 0,
        'exact_pos_close_seq': 0,
        'total': 0,                 # total based on reference entries with exact seq match
    }

    info_lines_printed = False
    lnum = 0
    
    with gzip.open(args.freeze_vcf, 'rt') as fh:
        cr = csv.reader(fh, delimiter='\t')
        for row in cr:
            lnum += 1
            if re.match(r'^\#.*$', row[0]):
                # take contig and FORMAT lines from freeze file
                m = re.match(r'^\#\#(contig|FORMAT)=.*$', row[0])
                if m or row[0] == '#CHROM':
                    # insert INFO lines from individual VCFs before FORMAT
                    if m and m.group(1) == 'FORMAT' and not info_lines_printed:
                        ofh.write("\n".join(info_lines) + "\n")
                        info_lines_printed = True
                    ofh.write("\t".join(row) + "\n")

                # echo some INFO lines from freeze file

                # INFO lines from freeze file:
                ##INFO=<ID=ID,Number=A,Type=String,Description="ID of merged variant set (max one variant per sample)">
                ##INFO=<ID=VARTYPE,Number=A,Type=String,Description="Variant class">
                ##INFO=<ID=SVTYPE,Number=A,Type=String,Description="Variant type">
                ##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length of ref and alt alleles">
                ##INFO=<ID=SAMPLE,Number=.,Type=String,Description="Lead sample variant was called on. SV sequence, breakpoints, and contig locations come from the SV call in this sample.">
                ##INFO=<ID=REF_SD,Number=.,Type=Float,Description="Max segmental duplication (SD) identity variant intersects.">
                ##INFO=<ID=REF_TRF,Number=0,Type=Flag,Description="Variant intersects a reference TRF region">
                
                # already have ID (as PAV_ID), also keeping SAMPLE
                
                if re.match(r'^\#\#INFO=\<ID=(SAMPLE),.*$', row[0]):
                    ofh.write("\t".join(row) + "\n")
                continue
                    
            (chrom, pos, vcf_id, ref, alt, qual, filt, inf, fmt, *rest) = row
            ins_seq = alt[1:]
            
            # examine insertions only
            if not re.match(r'.*-INS-.*', vcf_id):
                continue

            # parse INFO fields
            inf_d = {}
            for i in inf.split(';'):
                k = i
                v = None
                
                if re.match(r'^.*=.*$', i):
                    (k, v) = i.split('=')
                    
                if k in inf_d:
                    log_fatal("key " + k + " already seen in " + info)
                inf_d[k] = v
            
            # count 1s to find number of haplotype/samples matches
            n_haplotypes = 0
            for gt in rest:
                for hap in gt.split('|'):
                    if hap == '1':
                        n_haplotypes += 1
            
            # different levels of matching

            # exact: ref pos, insertion length, insertion sequence
            m_exact = []
            m_exact_d = {}
            
            # exact: insertion length, insertion sequence
            # close: ref pos
            m_exact_seq_close_pos = []
            m_exact_seq_close_pos_d = {}

            # exact: ref pos
            # close: insertion length, insertion sequence
            m_exact_pos_close_seq = []
            m_exact_pos_close_seq_d = {}

            # -----------------------------------------------------------
            # find matching MEIs by sequence position + ins length
            # -----------------------------------------------------------
            ref_pos_in_index = False
            ins_len = len(ins_seq)
            ref_pos = ":".join([chrom, pos, str(ins_len)])
            
            if ref_pos in rpi:
                ref_pos_in_index = True
                meis = rpi[ref_pos]
                
                # check for exact/approximate sequence match
                num_em = 0
                for mei in meis:
                    (m_chrom, m_pos, m_vcf_id, m_ref, m_alt, m_qual, m_flt, m_inf, m_fmt, *m_rest) = mei
                    m_ins_seq = m_alt[1:]

                    # case 1: exact position match and exact sequence match
                    if m_ins_seq == ins_seq:
                        num_em += 1
                        m_exact.append(mei)
                        m_exact_d[m_vcf_id] = True

                    # case 2: exact position match and close sequence match
                    # TODO
                    elif False:
                        pass

            # -----------------------------------------------------------
            # find matching MEIs by sequence identity
            # -----------------------------------------------------------
            ref_seq_in_index = False
            ins_seq = alt[1:]
            
            if ins_seq in si:
                ref_seq_in_index = True
                meis = si[ins_seq]

                # filter out the exact matches that have already been found
                new_meis = []
                for mei in meis:
                    (m_chrom, m_pos, m_vcf_id, m_ref, m_alt, m_qual, m_flt, m_inf, m_fmt, *m_rest) = mei
                    if m_vcf_id not in m_exact_d:
                        new_meis.append(mei)

                n_new = len(new_meis)

                # check exact sequence matches for proximity to ref
                if n_new > 0:
                    for new_mei in new_meis:
                        (m_chrom, m_pos, m_vcf_id, m_ref, m_alt, m_qual, m_flt, m_inf, m_fmt, *m_rest) = new_mei
                        if chrom != m_chrom:
                            continue
                        dist = abs(int(pos) - int(m_pos))
                        if dist < args.max_dist_bp:
                            m_exact_seq_close_pos.append(new_mei)
                            m_exact_seq_close_pos_d[m_vcf_id] = True
                        
                n_exact_seq_close_pos = len(m_exact_seq_close_pos)

            # -----------------------------------------------------------
            # output exemplar
            # -----------------------------------------------------------
            exact_seq_matches = m_exact
            exact_seq_matches.extend(m_exact_seq_close_pos)

            if len(exact_seq_matches) > 0:
                # index exact matches by sample-haplotype
                esm_d = {}
                for mei in exact_seq_matches:
                    (m_chrom, m_pos, m_vcf_id, m_ref, m_alt, m_qual, m_flt, m_inf, m_fmt, *m_rest) = mei
                    (sample, hap, m_id) = m_vcf_id.split(":")
                    key = "-".join([sample, hap])
                    if key in esm_d:
                        fatal("duplicate key " + key)
                    esm_d[key] = mei
                    
                # freeze SAMPLE (e.g., 'NA18534-h1' should specify the corresponding sample/haplotype)
                sample = inf_d['SAMPLE']
                ind_mei = esm_d[sample]
                (i_chrom, i_pos, i_vcf_id, i_ref, i_alt, i_qual, i_flt, i_inf, i_fmt, *i_rest) = ind_mei

                # sanity checks
                if i_ref != ref:
                    fatal("REF mismatch")
                if i_alt != alt:
                    fatal("ALT mismatch")

                (samp, hap, m_id) = i_vcf_id.split(":")
                i_inf += ";SAMPLE=" + sample
                
                # keep position from freeze VCF; individual MEI location may be different
                new_cols = [chrom, pos, m_id, ref, alt, qual, filt, i_inf, fmt, *rest]
                
                # also need haplotypes from freeze
                ofh.write("\t".join(new_cols) + "\n")

                # record type
                m = re.match(r'^.*-(ALU|SVA|LINE1)$', m_id)
                if not m:
                    fatal("unable to parse ME type from id " + m_id)
                ME_types_found[m.group(1)] += 1

            # -----------------------------------------------------------
            # update counts
            # -----------------------------------------------------------
            if ref_pos_in_index:
                n_freeze_vcf_pos_plus_inslen_found += 1
            if ref_seq_in_index:
                n_freeze_vcf_ins_seqs_found += 1
            if ref_pos_in_index or ref_seq_in_index:
                entry_counts['total'] += 1
                hap_counts['total'] += n_haplotypes;
                freeze_vcf_pos_found[chrom + ":" + pos] = True
                
            n_exact = len(m_exact)
            n_exact_seq_close_pos = len(m_exact_seq_close_pos)
            n_exact_pos_close_seq = len(m_exact_pos_close_seq)

            hap_counts['exact'] += n_exact;
            hap_counts['exact_seq_close_pos'] += n_exact_seq_close_pos;
            hap_counts['exact_pos_close_seq'] += n_exact_pos_close_seq;

            # entry counts are nonredundant
            if n_exact > 0:
                entry_counts['exact'] += 1;
            elif n_exact_seq_close_pos > 0:
                entry_counts['exact_seq_close_pos'] += 1;
            elif n_exact_pos_close_seq > 0:
                entry_counts['exact_pos_close_seq'] += 1;

    # -----------------------------------------------------------
    # report stats
    # -----------------------------------------------------------
    info("freeze VCF ref positions matched : " + str(len(freeze_vcf_pos_found)))
    info("freeze VCF lines matched by position and SVLEN : " + str(n_freeze_vcf_pos_plus_inslen_found))
    info("freeze VCF lines matched by insertion sequence : " + str(n_freeze_vcf_ins_seqs_found))

    for key in ['total', 'exact', 'exact_seq_close_pos', 'exact_pos_close_seq']:
        info("entries / " + key + " : " + str(entry_counts[key]))

    for key in ['total', 'exact', 'exact_seq_close_pos', 'exact_pos_close_seq']:
        info("haplotypes / " + key + " : " + str(hap_counts[key]))

    for key in ME_types_found.keys():
        info(key + " : " + str(ME_types_found[key]))
        
    # TODO - check for matches that are close in sequence and position but not identical in either?
    # TODO - check that all the exact sequence matches match (except for vcf_id, which contains sample + hap)
    
    # TODO - for each position add to the VCF:
    #   -actual number of exact occurrences (i.e., haplotype+sample count, with or without 'un')
    #   -number of merged occurrences (i.e., count up the '1's in the reported genotypes)

    # total hap count (all 1s)
    # exact hap count
    # other categories of counts
    
    # TODO - check that the accumulated matches are found in the expected samples/haplotypes [optional]
    
if __name__ == '__main__':
    main()
