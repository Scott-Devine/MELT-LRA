#!/usr/bin/env python3

# Initial version of MELT-RISC pipeline.

import argparse
import csv
import gzip
import hashlib
import os
import re
import sys

# ------------------------------------------------------
# Globals
# ------------------------------------------------------
FASTA_SUFFIX_RE = r'\.(fa.gz|fasta.gz)$'
FASTA_FILE_RE = r'^(.*)' + FASTA_SUFFIX_RE
DEBUG = False

# increase CSV max field size
csv.field_size_limit(256 * 1024 * 1024)

MIN_TSD_LEN = 1
MAX_TSD_LEN = 40

# perl -e 'while (<>) { chomp; if (!(/^>/)) { $l += length($_);} } print "length=$l\n"' <SVA_A.fa 
ME_LENGTHS = {
    'ALU': 281,
    'SVA': 1316,    # MELT SVA reference
    'SVA_A': 1387,
    'SVA_F': 1375,
    'LINE1': 6019
}

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
# read_single_fasta_file
# ------------------------------------------------------
def read_single_fasta_file(fpath):
    defline = None
    seq = ''
    with gzip.open(fpath, 'rt') as fh:
        for line in fh:
            # defline
            m = re.match(r'^>(.*)$', line)
            if m:
                if defline is not None:
                    fatal("multiple deflines found in single sequence FASTA file")
                defline = m.group(1)
            else:
                rsl = line.rstrip()
                seq = seq + rsl
                rsl = rsl.upper()
                rsl = re.sub(r'[\s\-]+', '', rsl)

    seqlen = len(seq)
    info("read sequence " + defline + " of length " + str(seqlen) + " from " + fpath)
    return { 'defline': defline, 'seq': seq, 'len': seqlen }

# ------------------------------------------------------
# read_fasta_dir
# ------------------------------------------------------
def read_fasta_dir(dpath, seqid, skip_seqids):
    files = {}

    for file in os.listdir(dpath):
        m = re.match(FASTA_FILE_RE, file)
        if m:
            # read specified reference sequence only
            if ((seqid is not None) and (m.group(1) != seqid)) or (m.group(1) in skip_seqids):
                continue
            file_id = m.group(1)
            fpath = os.path.join(dpath, file)
            ff = read_single_fasta_file(fpath)
            if file_id in files:
                fatal("duplicate FASTA file id " + file_id)
            files[file_id] = ff
    return files
        
# ------------------------------------------------------
# read_vcf_contigs
# ------------------------------------------------------
def read_vcf_contigs(vpath):
    contigs_l = []
    contigs_d = {}
    lnum = 0
    with gzip.open(vpath, 'rt') as fh:
        for line in fh:
            lnum += 1
            m = re.match(r'^\#\#contig=\<ID=([^,]+),length=(\d+),md5=([a-z0-9]{32})\>$', line)
            if m:
                c_id = m.group(1)
                c_len = int(m.group(2))
                c_md5 = m.group(3)
                debug("file id=" + c_id + " len=" + str(c_len) + " md5=" + c_md5)
                contig = { 'id': c_id, 'len': c_len, 'md5': c_md5 }
                contigs_l.append(contig)
                if id in contigs_d:
                    fatal("duplicate contig id " + c_id + " at line " + str(lnum))
                contigs_d[c_id] = contig

    n_contigs = len(contigs_l)
    msg ="read " + str(n_contigs) + " contig(s) from " + vpath 
    info(msg)
    return { 'contigs_l': contigs_l, 'contigs_d': contigs_d }

# ------------------------------------------------------
# read_vcf_insertions
# ------------------------------------------------------
def read_vcf_insertions(vpath, ref_seqs, seqid, skip_seqids):
    insertions = []
    lnum = 0

    with gzip.open(vpath, 'rt') as fh:
        cr = csv.reader(fh, delimiter='\t')
        for row in cr:
            if re.match(r'^#', row[0]):
                continue
            
            (chrom, pos, vcf_id, ref, alt, qual, filt, info, fmt, gt) = row
            pos = int(pos)

            if ((seqid is not None) and (chrom != seqid)) or (chrom in skip_seqids):
                continue
            
            inf_d = {}
            for inf in info.split(';'):
                (k, v) = inf.split('=')
                if k in inf_d:
                    fatal("key " + k + " already seen in " + info)
                inf_d[k] = v

            if inf_d['SVTYPE'] != 'INS':
                continue

            # sanity check - reference sequence is a single base
            if len(ref) != 1:
                fatal("len(REF) != 1")

            # sanity check - reference base equals first base of alt sequence
            if alt[0] != ref:
                fatal("ALT[0] (" + alt[0] + " != REF (" + ref + ")")

            # sanity check - check reference base against reference sequence
            refseq = ref_seqs[chrom]
            ref_base = refseq['seq'][pos-1]
            if ref_base.upper() != ref:
                fatal("reference file base at " + str(pos) + " (" + ref_base  + ") != REF (" + ref + ")")
                
            insertions.append({
                'chrom': chrom,
                'pos': pos,
                'vcf_id': vcf_id,
                'ref': ref,
                'alt': alt,
                'len': len(alt),
                'ins': alt[1:],
                'qual': qual,
                'filt': filt,
                'info': info,
                'fmt': fmt,
                'gt': gt
            })
            
    return insertions
                
# ------------------------------------------------------
# read_water()
# ------------------------------------------------------
def read_water(wfile):
    n_lines = 0
    indel_lengths = {}
    alignment = None
    # alignments indexed by insertion seq name e.g., chr1-191377-INS-2012
    alignments = {}

    def process_alignment(al):
        vn = al['variant_name']
        if vn in alignments:
            fatal("duplicate variant name " + vn)
        alignments[vn] = alignment

        def update_coords(coords, mline):
            m = re.match(r'^\S+\s+(\d+) \S+\s+(\d+)\s*$', mline)
            if not m:
                fatal("unable to parse match line '" + mline + "'")
            c1 = int(m.group(1))
            c2 = int(m.group(2))
            (l, h) = (c1, c2) if c1 <= c2 else (c2, c1)
            if coords['x1'] is None or coords['x1'] > l:
                coords['x1'] = l
            if coords['x2'] is None or coords['x2'] < h:
                coords['x2'] = h
        
        # TODO - parse min/max coords from match lines
        me_coords = { 'x1': None, 'x2': None }
        insertion_coords = { 'x1': None, 'x2': None }
        ctr = -1

        for ml in al['match_lines']:
            ctr = (ctr + 1) % 4
            if ctr == 1:
                if re.match(r'^\s*$', ml):
                    continue
                update_coords(me_coords, ml)
            elif ctr == 3:
                update_coords(insertion_coords, ml)
        
        al['ME_coords'] = me_coords
        al['insertion_coords'] = insertion_coords
        # save memory
        al['match_lines'] = None
        
    # reading detailed alignment match lines
    reading_match_lines = False
        
    with gzip.open(wfile, 'rt') as fh:
        for line in fh:
            n_lines += 1
            if re.match(r'^\# Aligned_sequences: 2.*$', line):
                if alignment is not None:
                    process_alignment(alignment)
                alignment = { 'match_lines': [] }
            else:
                m = re.match(r'^# 1: (\S+).*', line)
                if m:
                    alignment['ME'] = m.group(1)
                    
                m = re.match(r'^# 2: (\S+-(INS|DEL)-(\d+))', line)
                if m:
                    alignment['variant_name'] = m.group(1)
                    alignment['variant_type'] = m.group(2)
                    alignment['variant_len'] = int(m.group(3))
                    if alignment['variant_len'] < 50:
                        print("lnum=" + str(n_lines))

                m = re.match(r'^# Identity:\s+(\d+)\/(\d+) \((\s*[\d\.]+)%\)', line)
                if m:
                    alignment['matches'] = int(m.group(1))
                    alignment['length'] = int(m.group(2))
                    alignment['pct_id'] = float(m.group(3))

                m = re.match(r'^# Gaps:\s+(\d+)\/(\d+) \((\s*[\d\.]+)%\)', line)
                if m:
                    alignment['gaps'] = int(m.group(1))

                m = re.match(r'^# Score: ([\d\.]+).*$', line)
                if m:
                    alignment['score'] = float(m.group(1))

                m = re.match(r'^#=======================================\s*$', line)
                if alignment is not None and 'ME' in alignment and m:
                    reading_match_lines = True if not reading_match_lines else False
                elif re.match(r'^#\-+\s*$', line):
                    pass
                elif reading_match_lines:
                    alignment['match_lines'].append(line)
                    
    process_alignment(alignment)
    return alignments

# ------------------------------------------------------
# check_insertion_for_tsd
# ------------------------------------------------------
def check_insertion_for_tsd(ins, ref_seqs):
    ref_seq = ref_seqs[ins['chrom']]
    ref_seq_len = ref_seq['len']
    ins_seq = ins['ins']
    ins_seq_len = len(ins_seq)

    # find longest TSD _after_ the insertion

    # forward strand:
    # [<TSD>........polyA]<TSD>
    #
    # reverse strand:
    # [<TSD>polyT........]<TSD>
    #
    # [] = inserted sequence

    # TODO - note that this situation is possible too, but it depends on the insertion caller
    # where the 5' TSD appears wrt to the called insertion. PAV always seems to place the
    # 5' TSD inside the insertion
    #
    # <TSD>][........polyA<TSD>]
    
    ref_seq_pos = ins['pos']
    ins_seq_pos = 0
    tsd_after = ''

    while ref_seq['seq'][ref_seq_pos].upper() == ins_seq[ins_seq_pos].upper():
        tsd_after = tsd_after + ref_seq['seq'][ref_seq_pos]
        ref_seq_pos += 1
        ins_seq_pos += 1
        if ref_seq_pos >= ref_seq_len or ins_seq_pos >= ins_seq_len:
            break

    # find longest TSD _before_ the insertion:
    # <TSD>][........polyA<TSD>]
    ref_seq_pos = ins['pos'] - 1
    ins_seq_pos = ins_seq_len - 1
    tsd_before = ''
    while ref_seq['seq'][ref_seq_pos].upper() == ins_seq[ins_seq_pos].upper():
        tsd_before = ref_seq['seq'][ref_seq_pos] + tsd_before
        ref_seq_pos -= 1
        ins_seq_pos -= 1
        if ref_seq_pos < 0 or ins_seq_pos < 0:
            break

    return {
        'before': { 'tsd': tsd_before, 'len': len(tsd_before) },
        'after': { 'tsd': tsd_after, 'len': len(tsd_after) },
    }

def check_insertion_for_ME_match(vcf_ins, aligns, min_pctid, min_pctid_nogaps, min_pctcov, strand):
    ins_name = vcf_ins['chrom'] + '-' + str(vcf_ins['pos']+1) + '-INS-' + str(vcf_ins['len']-1)
    if ins_name not in aligns:
        # ugh
        ins_name = vcf_ins['chrom'] + '-' + str(vcf_ins['pos']) + '-INS-' + str(vcf_ins['len']-1)
        if ins_name not in aligns:
            fatal("no alignment found for " + ins_name)
    al = aligns[ins_name]
    me_match = None

    # compute percent identity without gaps
    pctid_nogaps = (al['matches'] / (al['length'] - al['gaps'])) * 100.0
    
    debug("checking " + ins_name + " for match, al=" + str(al))
    if (al['pct_id'] >= min_pctid) and (pctid_nogaps >= min_pctid_nogaps):
        debug("checking " + ins_name + " for match, pctid is good ")
        # TODO - pctcov shouldn't include TSD or polyA
#        available_bp = vcf_ins['len'] - vcf_ins['tsd']['len']
        available_bp = vcf_ins['len']
        ins_pctcov = (al['length'] / available_bp) * 100.0
        me_pctcov =  (al['length'] / ME_LENGTHS[al['ME']]) * 100.0
        debug("checking " + ins_name + " for match, ins_pctcov=" + str(ins_pctcov) + " me_pctcov=" + str(me_pctcov))
        if ins_pctcov >= min_pctcov:
            debug("checking " + ins_name + " for match, ins_pctcov is good, setting match to nonempty")
            me_match = { 'ME': al['ME'], 'pctid':al['pct_id'], 'ins_pctcov': ins_pctcov, 'me_pctcov': me_pctcov, 'alignment': al, 'strand': strand }
    return me_match

def check_insertion_for_polyA(vcf_ins):
    ins = vcf_ins['ins']
    
    # TODO - don't require exact match
    def find_polyX(base, start, offset):
        end = start
        plen = 0
        
        while (end >= 0) and (end < len(ins)) and ins[end] == base:
            end += offset
            plen += 1

        x1 = start if offset > 0 else end
        x2 = end if offset > 0 else start
        return { 'len': plen, 'x1': x1, 'x2': x2 }
    
    # 1. search backwards from end of insertion sequence for polyA
    vcf_ins['polyA'] = find_polyX('A', len(ins) - 1, -1)
    # 2. search forwards from start of insertion sequence for polyT
    vcf_ins['polyT'] = find_polyX('T', vcf_ins['tsds']['after']['len'], 1)

# ------------------------------------------------------
# print_insertion()
# ------------------------------------------------------
def print_insertion(ins, ref_seqs):
    ref_seq = ref_seqs[ins['chrom']]
    ref_seq_len = ref_seq['len']
    ref_seq_pos = ins['pos']
    l_context_bp = 10
    r_context_bp = 30
    
    # display region around insertion point (i.e., the point after the REF base)
    # TODO - check for off by 1 errors:
    bp_before = l_context_bp if ref_seq_pos > l_context_bp else ref_seq_pos
    bp_after = r_context_bp if (ref_seq_pos + r_context_bp < ref_seq_len) else ref_seq_len - ref_seq_pos

    sb_from = ref_seq_pos - bp_before
    seq_before = ref_seq['seq'][sb_from:ref_seq_pos]
    sa_to = ref_seq_pos + bp_after
    seq_after = ref_seq['seq'][ref_seq_pos:sa_to]

    # last base of seq_before should be REF
    if seq_before[-1].upper() != ins['ref']:
        fatal("seq_before[-1] (" + seq_before[-1] + ") != REF (" + ins['ref'] + ")")

    tsd_str = ins['tsds']['after']['tsd'] if ins['tsds']['after']['tsd'] is not None else '-'
    pos_str = (ins['chrom'] + ":" + str(ins['pos'])).ljust(16)
    # sequence to display inside the insertion
    l_bp = 45
    r_bp = 45
    middle_bp = 14
    remaining_bp = ins['len'] - 1 - l_bp - r_bp
    ins_str = ins['alt'][1:l_bp+1] + "..." + ("+" + str(remaining_bp) + "bp").center(middle_bp-6,".") + "..." + ins['alt'][-r_bp:]

    pctid_str = (str(ins['me_match']['pctid']) + "%").rjust(6)
    # percent of insert covered by alignment
    ins_pctcov_str = (("%.1f" % ins['me_match']['ins_pctcov']) + "%").rjust(6)
    # percent of reference ME sequence covered by alignment
    me_pctcov_str = (("%.1f" % ins['me_match']['me_pctcov']) + "%").rjust(6)
    me_str = ins['me_match']['ME'].ljust(5) + "|" + ins['me_match']['strand'].ljust(3) + "|" + me_pctcov_str + "|" + pctid_str + "|" + ins_pctcov_str
    print(pos_str + "|" + me_str + "| " + seq_before + " [" + ins_str + "] " + seq_after)

    # print TSD, polyA position
    pos_str = re.sub(r'\S', ' ', pos_str)
    me_str = re.sub(r'\S', ' ', me_str)
    seq_before = re.sub(r'\S', ' ', seq_before)
    seq_after = ""
    isl = len(ins_str)
    ins_str_left = "".center(l_bp)
    ins_str_right = "".center(r_bp)
    ins_str_middle = "".center(middle_bp)
    
    # TSD
    tl = ins['tsds']['after']['len']
    if tl > 0:
        rep = '^' + 'TSD'.center(tl-2, '^') + '^'
        rep = rep[0:tl]
        ins_str_left = rep.ljust(l_bp)[0:l_bp]
        seq_after = rep

    # polyA
    pa = ins['polyA']
    pa_len = pa['len']
    if pa_len >= 7:
        ins_str_right = ('<' + "polyA".center(pa_len-2, "-") + '>').rjust(r_bp)[-r_bp:]
        
    # polyT
    pt = ins['polyT']
    pt_len = pt['len']
    if pt_len >= 7:
        ins_str_left = ins_str_left[0:pt['x1']] + '<' + "polyT".center(pt_len-2, "-") + '>' + ins_str_left[pt['x2']:]
        ins_str_left = ins_str_left[0:l_bp]
        
    print(pos_str + " " + me_str + "  " + seq_before + "  " + ins_str_left + ins_str_middle + ins_str_right + "  " + seq_after)

    # print ME alignment
    ins_coords = ins['me_match']['alignment']['insertion_coords']
    x1 = ins_coords['x1']
    x2 = ins_coords['x2']
    ins_ld = "<" if ins['me_match']['strand'] == '-' else '['
    ins_rd = "]" if ins['me_match']['strand'] == '-' else '>'
    ins_len = len(ins['ins'])
    ins_left = "".center(x1-1) + ins_ld + ins['me_match']['ME'].ljust(x2 - x1 + 1, "-") + ins_rd if x1 < l_bp else "".center(l_bp)
    ins_right = ins_ld + ins['me_match']['ME'].rjust(x2 - x1 + 1, "-") + ins_rd + "".center(ins_len - x2) if x2 > (ins_len - r_bp) else "".center(r_bp)
    ins_str = ins_left[0:l_bp] + "".center(14) + ins_right[-r_bp:]
    print(pos_str + " " + me_str + "  " + seq_before + "  " + ins_str)
    # DEBUG
    if ins['tsds']['before']['tsd'] != "":
        print("TSD before=" + ins['tsds']['before']['tsd'])
    print()
    
# ------------------------------------------------------
# main()
# ------------------------------------------------------
def main():

    # input
    parser = argparse.ArgumentParser(description='Check FASTA sequences against a set of VCF reference contigs.')
    parser.add_argument('--vcf', required=True, help='Path to VCF file containing contigs to check.')
    parser.add_argument('--fasta_dir', required=True, help='Path to directory that contains FASTA reference files.')
    parser.add_argument('--alu_water', required=True, help='Path to ALU water alignment output file.')
    parser.add_argument('--alu_water_rev', required=True, help='Path to ALU water reverse strand alignment output file.')
    parser.add_argument('--sva_water', required=True, help='Path to SVA water alignment output file.')
    parser.add_argument('--sva_water_rev', required=True, help='Path to SVA water reverse strand alignment output file.')
    parser.add_argument('--line_water', required=True, help='Path to LINE1 water alignment output file.')
    parser.add_argument('--line_water_rev', required=True, help='Path to LINE1 water reverse strand alignment output file.')
    parser.add_argument('--min_seqlen', required=False, type=int, default=100, help='Minimum insertion sequence length.')
    parser.add_argument('--min_pctid', required=False, type=int, default=90, help='Minimum percent identity of alignment.')
    parser.add_argument('--min_pctid_nogaps', required=False, type=int, default=90, help='Minimum percent identity of alignment ignoring gaps.')
    parser.add_argument('--min_pctcov', required=False, type=int, default=85, help='Minimum percent coverage of alignment.')
    parser.add_argument('--seqid', required=False, help='Optional sequence id: process only insertions on this reference sequence.')
    parser.add_argument('--skip_seqids', required=False, help='Optional comma-delimited list of sequence ids to skip.')
    args = parser.parse_args()

    skip_seqids = {}
    if args.skip_seqids is not None and args.skip_seqids != '':
        for seqid in args.skip_seqids.split(','):
            info("skipping insertions on sequence " + seqid)
            skip_seqids[seqid] = True
    
    # read reference FASTA files
    fasta_files = read_fasta_dir(args.fasta_dir, args.seqid, skip_seqids)

    # read VCF contigs
    vcf_contigs = read_vcf_contigs(args.vcf)

    # TODO - incorporate MD5 check from checkVCFReferenceSeqs.py
    
    # read VCF insertions
    vcf_insertions = read_vcf_insertions(args.vcf, fasta_files, args.seqid, skip_seqids)
    info("read " + str(len(vcf_insertions)) + " insertions from " + args.vcf)
    
    # filter by length
    vcf_insertions = [i for i in vcf_insertions if i['len'] >= args.min_seqlen]
    info("read " + str(len(vcf_insertions)) + " insertions of length >= " + str(args.min_seqlen))

    # read alignment files
    alu_aligns = read_water(args.alu_water)
    info("read " + str(len(alu_aligns)) + " ALU alignment(s) from " + args.alu_water)
    alu_rev_aligns = read_water(args.alu_water_rev)
    info("read " + str(len(alu_rev_aligns)) + " reverse strand ALU alignment(s) from " + args.alu_water_rev)

    sva_aligns = read_water(args.sva_water)
    info("read " + str(len(sva_aligns)) + " SVA alignment(s) from " + args.sva_water)
    sva_rev_aligns = read_water(args.sva_water_rev)
    info("read " + str(len(sva_rev_aligns)) + " reverse strand SVA alignment(s) from " + args.sva_water_rev)

    line_aligns = read_water(args.line_water)
    info("read " + str(len(line_aligns)) + " LINE1 alignment(s) from " + args.line_water)
    line_rev_aligns = read_water(args.line_water_rev)
    info("read " + str(len(line_rev_aligns)) + " reverse strand LINE1 alignment(s) from " + args.line_water_rev)

    print("Location        |ME   |+/-|%ME   |%id   |%cov  | insertion")

    n_tsd_before = 0
    n_tsd_after = 0
    n_me_match = 0
    for vcf_ins in vcf_insertions:
        # check for presence of target site duplication
        tsds = check_insertion_for_tsd(vcf_ins, fasta_files)
        vcf_ins['tsds'] = tsds
        
        if tsds['before']['len'] >= MIN_TSD_LEN and tsds['before']['len'] <= MAX_TSD_LEN:
            n_tsd_before += 1
        if tsds['after']['len'] >= MIN_TSD_LEN and tsds['after']['len'] <= MAX_TSD_LEN:
            n_tsd_after += 1

        # check for qualifying match with mobile element
        alu_match = check_insertion_for_ME_match(vcf_ins, alu_aligns, args.min_pctid, args.min_pctid_nogaps, args.min_pctcov, '+')
        alu_rev_match = check_insertion_for_ME_match(vcf_ins, alu_rev_aligns, args.min_pctid, args.min_pctid_nogaps, args.min_pctcov, '-')
        sva_match = check_insertion_for_ME_match(vcf_ins, sva_aligns, args.min_pctid, args.min_pctid_nogaps, args.min_pctcov, '+')
        sva_rev_match = check_insertion_for_ME_match(vcf_ins, sva_rev_aligns, args.min_pctid, args.min_pctid_nogaps, args.min_pctcov, '-')
        line_match = check_insertion_for_ME_match(vcf_ins, line_aligns, args.min_pctid, args.min_pctid_nogaps, args.min_pctcov, '+')
        line_rev_match = check_insertion_for_ME_match(vcf_ins, line_rev_aligns, args.min_pctid, args.min_pctid_nogaps, args.min_pctcov, '-')
        matches = [m for m in [alu_match, alu_rev_match, sva_match, sva_rev_match, line_match, line_rev_match] if m is not None]
        n_matches = len(matches)

        # sort matches and pick the highest-scoring
        me_match = None
        if n_matches > 0:
            sorted_matches = sorted(matches, key = lambda x: x['alignment']['score'], reverse=True)
            me_match = sorted_matches[0]
            # kick the can down the road
            if n_matches > 1 and sorted_matches[0]['alignment']['score'] == sorted_matches[1]['alignment']['score']:
                fatal("multiple qualifying matches with the same score: " + str(sorted_matches))
                
        vcf_ins['me_match'] = me_match
        if me_match is not None:
            n_me_match += 1
            # don't expect this to happen for PAV-called L1-mediated insertions:
            if tsds['before']['len'] > tsds['after']['len']:
                fatal("Longer TSD sequence found _before_ the insertion point.")

        # check for polyA/polyT
        polyA = check_insertion_for_polyA(vcf_ins)

        # print insertions with ME match (but maybe no TSD)
        if me_match is not None:
            print_insertion(vcf_ins, fasta_files)
        
    # summary
    info("read " + str(len(vcf_insertions)) + " insertions from " + args.vcf)
    info("found " + str(n_tsd_before) + " insertion(s) of length >= " + str(args.min_seqlen) + " with TSDs _before_ the insertion point")
    info("found " + str(n_tsd_after) + " insertion(s) of length >= " + str(args.min_seqlen) + " with TSDs _after_ the insertion point")
    info("found " + str(n_me_match) + " insertion(s) of length >= " + str(args.min_seqlen) + " with ME matches")
    
if __name__ == '__main__':
    main()

