#!/usr/bin/env python3

# Convert MELT-RISC CSV output to VCF.

import argparse
import csv
from datetime import date
import gzip
import hashlib
import os
from pathlib import Path
import re
import sys

# ------------------------------------------------------
# Globals
# ------------------------------------------------------
VERSION = '1.4.2'

CSV_FILE_RE = r'^(.*)\.csv'
CSV_SAMPLE_ID_RE = r'^([^\/]+)-PAV-MEs.*$'
CSV_HEADER = None

DEBUG = False

# increase CSV max field size
csv.field_size_limit(256 * 1024 * 1024)

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
# read_csv()
# ------------------------------------------------------
def read_csv(csv_path):
    # list of MEIs in order encountered
    meis = []
    # map ref position to MEI
    pos2mei = {}
    
    with open(csv_path, 'r') as csv_fh:
        cr = csv.reader(csv_fh, delimiter=',')

        for row in cr:
            (samples, chrom, pos, strand, *rest) = row
            # skip header
            if samples == 'samples':
                continue
            
            key = ":".join([chrom, pos])
            if key in pos2mei:
                fatal()

            mei = { 'csv': row }
            meis.append(mei)
            pos2mei[key] = mei

    return { 'meis': meis, 'pos2mei': pos2mei }
            
# ------------------------------------------------------
# read_input_vcf()
# ------------------------------------------------------
def read_input_vcf(vcf_path, meis):
    pos2mei = meis['pos2mei']
    contigs = []
    with gzip.open(vcf_path, 'rt') as fh:
        for line in fh:

            # comment line
            if re.match(r'^##.*$', line):
                ## contig
                m = re.match(r'^##contig=<ID=([^,]+),length=(\d+),md5=([^>]+)>', line)
                if m:
                    contig = { 'id': m.group(1), 'length': int(m.group(2)), 'md5': m.group(3), 'vcf': line}
                    contigs.append(contig)

            # non-comment line
            else:
                ## INS variant e.g., 
                # chr1	94824	chr1-94825-INS-3	C	CTTT	.	PASS	ID=chr1-94825-INS-3;SVTYPE=INS;SVLEN=3;HAP=h1;COV_MEAN=1.0;COV_PROP=1.0;QRY_REGION=haplotype1-0000013:67048705-67048707;QRY_STRAND=-;CALL_SOURCE=CIGAR	GT	1
                m = re.match(r'', line)
                if m:
                    spl = line.split("\t")
                    (chrom, pos, id, ref, alt, qual, filter, inf, fmt, sample) = tuple(line.split("\t"))
                
                    # check if it matches an MEI in the CSV file
                    key = ":".join([chrom, pos])
                    if key in pos2mei:
                        vcf = {
                            'chrom': chrom,
                            'pos': pos,
                            'id': id,
                            'ref': ref,
                            'alt': alt,
                            'qual': qual,
                            'filter' : filter,
                            'info': inf,
                            'format': fmt,
                            'sample' : sample.rstrip()
                        }
                        pos2mei[key]['vcf'] = vcf

    # check that every MEI has been matched with the corresponding variant
    for mei in meis['meis']:
        if 'vcf' not in mei:
            fatal("no VCF insertion found for MEI " + str(mei['csv']))
                
    return { 'contigs': contigs }

# strip % sign, not permitted in VCF Float values
def strip_pct(s):
    m = re.match(r'^(.*)%$', s)
    if m:
        s = m.group(1)
    return s
    
# ------------------------------------------------------
# main()
# ------------------------------------------------------
def main():

    # input
    parser = argparse.ArgumentParser(description='Convert MELT-RISC CSV output to VCF.')
    parser.add_argument('--csv', required=True, help='Path to input MELT-RISC CSV file.')
    parser.add_argument('--vcf_input', required=True, help='Path to original PAV variant call VCF file.')
    parser.add_argument('--vcf_output', required=True, help='Path to output VCF file.')
    args = parser.parse_args()

    # extract sample id from csv filename
    csv_basename = os.path.basename(args.csv)
    m = re.match(CSV_SAMPLE_ID_RE, os.path.basename(csv_basename))
    if not m:
        fatal("couldn't parse sample id from " + csv_basename)
    sample = m.group(1)
    info("csv=" + csv_basename + " sample=" + sample)

    # read CSV
    meis = read_csv(args.csv)
    
    # read input VCF
    vcf = read_input_vcf(args.vcf_input, meis)
    
    with open(args.vcf_output, 'w') as vcf_fh:
        vcf_fh.write("##fileformat=VCFv4.2\n")
        vcf_fh.write("##fileDate=" + date.today().strftime("%Y%m%d") + "\n")
        vcf_fh.write("##source=MELT-RISC " + VERSION + "\n")
        vcf_fh.write("##reference=" + args.vcf_input + "\n")

        # write all input contigs, not just those with PAV calls
        ##contig=
        for contig in vcf['contigs']:
            vcf_fh.write(contig['vcf'])

        ##INFO=
        heading_info = [
            { 'id': "ID", 'num': "1", 'type': "String", 'descr': 'MEI ID' },
            { 'id': "PAV_ID", 'num': "1", 'type': "String", 'descr': 'ID of corresponding PAV insertion' },
            { 'id': "ME", 'num': "1", 'type': "String", 'descr': 'Mobile element type: ALU, LINE1, or SVA.' },
            { 'id': "PCT_ME", 'num': "1", 'type': "Float", 'descr': 'Percentage of the reference ME covered by S-W alignments.' },
            { 'id': "PCT_ID", 'num': "1", 'type': "Float", 'descr': 'Average percent identity of alignment(s) with the reference ME.' },
            { 'id': "PCT_COV", 'num': "1", 'type': "Float", 'descr': 'Percentage of the insertion sequence covered by S-W alignments to the reference ME, after subtracting TSD and polyA/polyT.' },
            { 'id': "TSD_SEQ", 'num': "1", 'type': "String", 'descr': 'Sequence of the target site duplication (TSD).' },
            { 'id': "POLYX_COORDS", 'num': "1", 'type': "String", 'descr': 'Coordinates of polyA/polyT sequence within the insertion in 1-based coordinates.' },
            { 'id': "ME_COORDS", 'num': "1", 'type': "String", 'descr': 'Coordinates of the ME within the insertion in 1-based coordinates.' },
            { 'id': "ME_MATCH", 'num': "1", 'type': "String", 'descr': 'Match string for Smith-Waterman alignment of insertion sequence against reference ME.' },
            { 'id': "ME_FAMILY", 'num': "1", 'type': "String", 'descr': 'CAlu/LINEU-called mobile element family, or SVA for SVAs.' },
            { 'id': "ME_SUBFAMILY", 'num': "1", 'type': "String", 'descr': 'CAlu/LINEU-called mobile element subfamily, or SVA for SVAs.' },
            { 'id': "ME_START", 'num': "1", 'type': "Integer", 'descr': 'CAlu/LINEU-determined start of the mobile element within the insertion. Empty for SVAs.' },
            { 'id': "ME_STOP", 'num': "1", 'type': "Integer", 'descr': 'CAlu/LINEU-determined end of the mobile element within the insertion. Empty for SVAs.'  },
            { 'id': "ME_NUM_DIAG_MATCHES", 'num': "1", 'type': "Integer", 'descr': 'CAlu/LINEU number of diagonal matches.' },
            { 'id': "ME_NUM_DIFFS", 'num': "1", 'type': "Integer", 'descr': 'CAlu/LINEU number of differences with respect to the reference ME.' },
            { 'id': "ME_DIFFS", 'num': "1", 'type': "String", 'descr': 'CAlu/LINEU report of differents with the reference ME.' },
            { 'id': "ME_OVERLAPPING_REPEATS", 'num': ".", 'type': "String", 'descr': 'Comma-delimited list of annotated RepeatMasker repeats that overlap the PAV insertion site.' },
            
            # genotype
            { 'column': 'FORMAT', 'id': "GT", 'num': "1", 'type': "String", 'descr': 'Genotype' },
        ]

        # TO ADD (back):
        #  strand
        #  insertion length
        
        for hi in heading_info:
            column = 'INFO'
            if 'column' in hi:
                column = hi['column']
            vcf_fh.write("##" + column + "=<ID=" + hi['id'] + ",Number=" + hi['num'] + ",Type=" + hi['type'] + ',Description="' + hi['descr'] + '"')
            vcf_fh.write(",Source=" + '"MELT-RISC"' + ",Version=" + '"' + VERSION + '">' + "\n")
        
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Variant length">
##INFO=<ID=HAP,Number=1,Type=String,Description="List of haplotype names variant was identified in">
##INFO=<ID=COV_MEAN,Number=.,Type=String,Description="Mean coverage for each haplotype under the whole variant (not just breakpoints, INFO/HAP order)">
##INFO=<ID=COV_PROP,Number=.,Type=String,Description="Proportion of reference bases under the whole variant with at least one aligned query (assembly sequence) (INFO/HAP order)">
##INFO=<ID=QRY_REGION,Number=.,Type=String,Description="Region of the query (assembly sequence) where this variant was found (chrom:pos-end, 1-based, closed coordinates, not BED) (INFO/HAP order)">
##INFO=<ID=QRY_STRAND,Number=.,Type=String,Description="Orientation of the aligned query (assembly sequence) at this site (INFO/HAP order)">
##INFO=<ID=CALL_SOURCE,Number=.,Type=String,Description="How variant was called - CIGAR, ALIGN_TRUNC, etc (INFO/HAP order)">
##INFO=<ID=INNER_REF,Number=.,Type=String,Description="Inversion inner breakpoint in reference coordinates (INFO/HAP order)">
##INFO=<ID=INNER_TIG,Number=.,Type=String,Description="Inversion inner breakpoint in contig coordinates (INFO/HAP order)">
##INFO=<ID=SEQ,Number=.,Type=String,Description="SV or indel sequence">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##ALT=<ID=INV,Description="Inversion">
##FILTER=<ID=PASS,Description="Variant passed filters">
##FILTER=<ID=QRY_FILTER,Description="Query filter region">
##FILTER=<ID=COMPOUND,Description="Inside larger variant">
##FILTER=<ID=SVLEN,Description="Variant size out of bounds">
##FILTER=<ID=TRIM,Description="Alignment trimming removed variant region">
        
        col_headings = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample]
        vcf_fh.write("#" + "\t".join(col_headings) + "\n")

        # write MEI calls to VCF
        for mei in meis['meis']:
            vcf = mei['vcf']

            (samples, chrom, pos, strand, ME, pct_ME, pct_id, pct_cov, iseq,
             left_flank_seq, right_flank_seq, TSD_seq,
             polyX_coords, ME_coords, insertion_coords, match_string,
             ME_family, ME_subfamily, ME_start, ME_stop, ME_num_diag_matches, ME_num_diffs, ME_diffs,
             overlapping_annots, genotype, hap1_region, hap2_region) = mei['csv']
            
            # Q/C - check that REF + ins seq = ALT
            alt_seq = vcf['ref'] + iseq
            if alt_seq != vcf['alt']:
                fatal("computed ALT (" + alt_seq + ") does not match original VCF ALT (" + vcf['alt'] + ") for " + vcf['id'])

            mei_id = vcf['id'] + "-" + ME

            # TODO - add INFO from CSV
            # TODO - add any standard INFO fields from VCF?
            
            # INFO fields
            infs = [
                ['ID', mei_id],
                ['PAV_ID', vcf['id']],
                ['ME', ME],
                ['PCT_ME', strip_pct(pct_ME)],
                ['PCT_ID', strip_pct(pct_id)],
                ['PCT_COV', strip_pct(pct_cov)],
                ['TSD_SEQ', TSD_seq],
                ['POLYX_COORDS', polyX_coords],
                ['ME_COORDS', ME_coords],
                ['ME_MATCH', match_string],
                ['ME_FAMILY', ME_family],
                ['ME_SUBFAMILY', ME_subfamily],
                ['ME_START', ME_start],
                ['ME_STOP', ME_stop],
                ['ME_NUM_DIAG_MATCHES', ME_num_diag_matches],
                ['ME_NUM_DIFFS', ME_num_diffs],
                ['ME_DIFFS', ME_diffs],
                ['ME_OVERLAPPING_REPEATS', overlapping_annots],
            ]
            inf_str = ";".join([inf[0] + "=" + ('.' if inf[1] == '' else inf[1]) for inf in infs])
            
            vcf_cols = [
                vcf['chrom'],   # CHROM
                vcf['pos'],     # POS
                mei_id,         # ID
                vcf['ref'],     # REF
                vcf['alt'],     # ALT
                '.',            # QUAL - VCF missing value
                'PASS',         # FILTER - assumed PASS
                inf_str,        # INFO
                'GT',           # FORMAT
                '1'             # SAMPLE
            ]
               
            vcf_fh.write("\t".join(vcf_cols))
            vcf_fh.write("\n")

                
if __name__ == '__main__':
    main()
    
