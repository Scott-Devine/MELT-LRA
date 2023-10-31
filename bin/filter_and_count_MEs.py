#!/usr/bin/env python3

# Filter and summarize output from multiple runs of find_MEs.py

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
CSV_FILE_RE = r'^(.*)\.csv'
CSV_SAMPLE_ID_RE = r'^(.*)-PAV-MEs.*$'
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
# read_csv_dir
# ------------------------------------------------------
def read_csv_dir(dpath):
    files = []
    for file in os.listdir(dpath):
        m = re.match(CSV_FILE_RE, file)
        if m:
            files.append(file)
    return files

# ------------------------------------------------------
# filter_and_index_csv_file
# ------------------------------------------------------
def filter_and_index_csv_file(csv_dir, output_dir, output_suffix, cfile, filters, unique_seq, unique_calu):
    global CSV_HEADER
    
    # sample index
    ind = {}
    
    # parse sample_id from CSV filename
    m = re.match(CSV_SAMPLE_ID_RE, cfile)
    if not m:
        fatal("couldn't parse sample_id from " + cfile + " with regex '" + CSV_SAMPLE_ID_RE + "'")
    sample_id = m.group(1)
        
    ipath = os.path.join(csv_dir, cfile)
    # construct output filename
    ofile = re.sub(r'\.csv$', '-' + output_suffix + '.csv', cfile)
    opath = os.path.join(output_dir, ofile)

    # check opath doesn't exist
    if os.path.exists(opath):
        fatal("output file " + opath + " already exists: please remove it and try again")
    
    # read from ipath, write to opath
    lnum = 0
    n_read = 0
    n_written = 0
    with open(ipath, 'rt') as ifh:
        cr = csv.reader(ifh, delimiter=',')

        with open(opath, 'wt') as ofh:
            for row in cr:
                lnum += 1

                (chrom, pos, strand, ME, pct_ME, pct_id, pct_cov, iseq,
                 left_flank_seq, right_flank_seq, TSD_seq,
                 polyX_coords, ME_coords, insertion_coords, match_string,
                 ME_family, ME_subfamily, ME_start, ME_stop, ME_num_diag_matches, ME_num_diffs, ME_diffs,
                 overlapping_annots, genotype, hap1_region, hap2_region) = row
                
                # header line
                if chrom == 'chrom':
                    if CSV_HEADER is None:
                        CSV_HEADER = ['samples']
                        CSV_HEADER.extend(row)
                    ofh.write(",".join(row) + "\n")
                    continue

                # MEI line
                n_read += 1

                # -------------------------------------------------
                # apply ME-specific filters
                # -------------------------------------------------
                filter = False
                me_filters = filters[ME]

                # -------------------------------------------------
                # filter by overlapping repeat_type
                # -------------------------------------------------
                # filter out anything with an overlapping_annot in this set:
                rtypes = me_filters['repeat_types']
                
                oas = [] if overlapping_annots == '' else overlapping_annots.split("|")

                for oa in oas:
                    debug("rtypes=" + str(rtypes))
                    debug("oa=" + str(oa))

                    # e.g., chr1:12677496:12677768:-:MLT1F2:LTR:ERVL-MaLR:0
                    (oa_chrom, oa_cstart, oa_cend, oa_strand, oa_rep_name, oa_rep_class, oa_rep_fam, oa_rep_left) = oa.split(':')
                    if oa_rep_fam in rtypes:
                        filter = True

                # min polyA/polyT length filter
                min_polya = me_filters['min_polya']

                pxc = [int(x) for x in polyX_coords.split('-')]
                polyx_bp = pxc[1] - pxc[0] + 1 if (pxc[1] - pxc[0]) >= 0 else 0
                
                if polyx_bp <  min_polya:
                    filter = True
                        
                # -------------------------------------------------
                # write MEI to output file and add to index
                # -------------------------------------------------
                if not filter:
                    ofh.write(",".join(row) + "\n")
                    n_written += 1
                    # index by location
                    kc = [chrom, pos, strand]
                    # and insertion sequence, if requested
                    if unique_seq:
                        kc.append(iseq)
                    # and cALU/LINEu calls
                    if unique_calu:
                        kc.append(ME_subfamily + ":" + ME_diffs)
                    key = ":".join(kc)
                    ind[key] = row
                    
    # print number of records in and out
    pct_written_str = "%.1f" % ((n_written/n_read) * 100.0)
    info("in:" + ipath + "  out:" + opath + "  wrote: " + str(n_written) + " / " + str(n_read) + " (" + pct_written_str + "%)")

    return(sample_id, ind)

# ------------------------------------------------------
# count_MEs
# ------------------------------------------------------
def count_MEs(me_ind):
    type_counts = {}
    fam_counts = {}
    subfam_counts = {}

    def count(cts, key):
        if key not in cts:
            cts[key] = 0
        cts[key] += 1
    
    for key in me_ind:
        me = me_ind[key]
        (chrom, pos, strand, ME, pct_ME, pct_id, pct_cov, iseq,
         left_flank_seq, right_flank_seq, TSD_seq,
         polyX_coords, ME_coords, insertion_coords, match_string,
         ME_family, ME_subfamily, ME_start, ME_stop, ME_num_diag_matches, ME_num_diffs, ME_diffs,
         overlapping_annots, genotype, hap1_region, hap2_region) = me
        
        count(type_counts, ME)
        count(type_counts, 'total')

        count(fam_counts, ME_family)
        count(fam_counts, 'total')

        count(subfam_counts, ME_subfamily)
        count(subfam_counts, 'total')

    return { 
        'ME_type' : type_counts,
        'ME_family': fam_counts,
        'ME_subfamily': subfam_counts
    }

# ------------------------------------------------------
# find_unique_MEs
# ------------------------------------------------------
def find_unique_MEs(fh, ME_inds, output_dir):
    sample_ids = sorted(ME_inds.keys())
    
    # map each ME to the samples in which it appears
    unique_MEs = {}
    
    for sid in sample_ids:
        s_ind = ME_inds[sid]
        for key in s_ind:
            if key not in unique_MEs:
                unique_MEs[key] = []
            unique_MEs[key].append(sid)

    # write sample count histogram to counts file
    sc_hist = {}
    for k in unique_MEs:
        sc = len(unique_MEs[k])
        if sc not in sc_hist:
            sc_hist[sc] = 0
        sc_hist[sc] += 1

    fh.write("\t".join(["num_samples", "num_MEIs"]) + "\n")
    keys = [k for k in sc_hist.keys()]
    for k in sorted(keys, key=lambda x: int(x), reverse=True):
        fh.write("\t".join([str(k), str(sc_hist[k])]) + "\n")
    fh.write("\t".join(['total', str(len(unique_MEs))]) + "\n")

    # generate MEI CSV files
    sample_sig_to_meis = {}
    for k in unique_MEs:
        me_sample_ids = sorted(unique_MEs[k])
        ssig = "_".join(me_sample_ids)
        if ssig not in sample_sig_to_meis:
            sample_sig_to_meis[ssig] = []
        sample_sig_to_meis[ssig].append(k)

    mei_files = {}
        
    for ssig in sample_sig_to_meis:
        mei_keys = sample_sig_to_meis[ssig]
        me_sample_ids = ssig.split("_")
        n_samples = len(me_sample_ids)
        ofile = str(n_samples) + "-samples.csv"
        write_mode = "wt"
        
        if ofile in mei_files:
            write_mode = "at"
        else:
            mei_files[ofile] = True
            
        opath = os.path.join(output_dir, ofile)
        with open(opath, write_mode) as ofh:
            # read header only once
            if write_mode == "wt":
                ofh.write(",".join(CSV_HEADER) + "\n")
            
            # loop over ref genome loci
            for k in mei_keys:
                # group samples at this locus by sequence then genotype
                sg_groups = {}
                for sid in me_sample_ids:
                    s_ind = ME_inds[sid]
                    mei = s_ind[k]
                    (chrom, pos, strand, ME, pct_ME, pct_id, pct_cov, iseq,
                     left_flank_seq, right_flank_seq, TSD_seq,
                     polyX_coords, ME_coords, insertion_coords, match_string,
                     ME_family, ME_subfamily, ME_start, ME_stop, ME_num_diag_matches, ME_num_diffs, ME_diffs,
                     overlapping_annots, genotype, hap1_region, hap2_region) = mei

                    if iseq not in sg_groups:
                        sg_groups[iseq] = { 'gt' : {}, 'n_samples' : 0 }
                    if genotype not in sg_groups[iseq]['gt']:
                        sg_groups[iseq]['gt'][genotype] = { 'meis': [], 'samples': [], 'n_samples': 0, 'seq': iseq, 'genotype': genotype }

                    sg_groups[iseq]['n_samples'] += 1
                    sg_groups[iseq]['gt'][genotype]['samples'].append(sid)
                    sg_groups[iseq]['gt'][genotype]['n_samples'] += 1
                    sg_groups[iseq]['gt'][genotype]['meis'].append(mei)

                # print sample groups
                for seq in sorted(sg_groups.keys(), key=lambda x: sg_groups[x]['n_samples'], reverse = True):
                    sample_str = ""

                    gtype_list = sorted(sg_groups[seq]['gt'].keys(), key = lambda x: sg_groups[seq]['gt'][x]['n_samples'], reverse = True)
                    for gtype in gtype_list:
                        v = sg_groups[seq]['gt'][gtype]
                        samples = v['samples']
                        sample_str += " ".join(samples) + " [" + gtype + "] "

                    sample_str += " - " + str(sg_groups[seq]['n_samples']) + " sample(s)"

                    # each unique sequence corresponds to an output row
                    row = [c for c in v['meis'][0]]
                    # set genotypes to a list, clear haplotype 1/2 info
                    row[-1] = ''
                    row[-2] = ''
                    row[-3] = 'multiple' if len(gtype_list) > 1 else gtype_list[0]
                    
                    ofh.write(sample_str + "," + ",".join(row) + "\n")  

def list_to_dict(l):
    d = {}
    if l is not None:
        for item in l.split(","):
            d[item] = True
    return d    

# ------------------------------------------------------
# main()
# ------------------------------------------------------
def main():

    # input
    parser = argparse.ArgumentParser(description='Filter and summarize output from multiple runs of find_MEs.py')
    parser.add_argument('--csv_dir', required=True, help='Path to directory containing CSV files produced by find_MEs.py.')
    parser.add_argument('--output_dir', required=True, help='Path to directory where filtered output files should be written.')
    parser.add_argument('--output_suffix', required=False, default='filtered', help='Suffix to append to output files.')
    parser.add_argument('--require_unique_sequence', required=False, action=argparse.BooleanOptionalAction, help='Whether to require unique sequence when merging MEIs.')
    parser.add_argument('--require_unique_calu_lineu', required=False, action=argparse.BooleanOptionalAction, help='Whether to require unique cALU/LINEu calls when merging MEIs.')
    # SVA filters
    parser.add_argument('--sva_excluded_repeat_types', required=False, help='Exclude/filter SVAs whose overlapping repeat type is in this list.')
    parser.add_argument('--sva_min_polyA_bp', required=False, default=0, help='Minimum SVA polyA/polyT length.')
    # Alu filters
    parser.add_argument('--alu_excluded_repeat_types', required=False, help='Exclude/filter Alus whose overlapping repeat type is in this list.')
    parser.add_argument('--alu_min_polyA_bp', required=False, default=0, help='Minimum Alu polyA/polyT length.')
    # LINE filters
    parser.add_argument('--line_excluded_repeat_types', required=False, help='Exclude/filter LINEs whose overlapping repeat type is in this list.')
    parser.add_argument('--line_min_polyA_bp', required=False, default=0, help='Minimum LINE polyA/polyT length.')
    # global filters
    # TODO
    args = parser.parse_args()

    # ------------------------------------------------------
    # filters
    # ------------------------------------------------------
    # excluded repeat types
    ex_sva_rep_types_d = list_to_dict(args.sva_excluded_repeat_types)
    ex_alu_rep_types_d = list_to_dict(args.alu_excluded_repeat_types)
    ex_line_rep_types_d = list_to_dict(args.line_excluded_repeat_types)

    # min polyA length
    sva_min_polya = int(args.sva_min_polyA_bp)
    alu_min_polya = int(args.alu_min_polyA_bp)
    line_min_polya = int(args.line_min_polyA_bp)
    
    filters = {
        'SVA': { 'repeat_types': ex_sva_rep_types_d, 'min_polya': sva_min_polya },
        'ALU' : { 'repeat_types': ex_alu_rep_types_d, 'min_polya': alu_min_polya },
        'LINE1': { 'repeat_types': ex_line_rep_types_d, 'min_polya': line_min_polya }
    }
        
    # read, filter, and index csv_files
    csv_files = read_csv_dir(args.csv_dir)
    sample_inds = {}
    for cf in csv_files:
        unique_seq = args.require_unique_sequence
        unique_calu = args.require_unique_calu_lineu
        (sample_id, index) = filter_and_index_csv_file(args.csv_dir, args.output_dir, args.output_suffix, cf, filters, unique_seq, unique_calu)
        sample_inds[sample_id] = index

    # ------------------------------------------------------
    # write summary/counts file
    # ------------------------------------------------------
    counts_file = "summary-counts.tsv"
    cpath = os.path.join(args.output_dir, counts_file)

    # write per-sample counts
    sample_ids = sorted(sample_inds.keys())
    sample_counts = {}

    for sid in sample_ids:
        sample_counts[sid] = count_MEs(sample_inds[sid])

    with open(cpath, "wt") as cfh:
        # different count types
        for type in ('ME_type', 'ME_family', 'ME_subfamily'):
            # column headings
            counts_headers = [type]
            counts_headers.extend(sample_ids)
            cfh.write("\t".join(counts_headers) + "\n")

            # take union of keys over all samples
            all_keys = {}
            for sid in sample_ids:
                for k in sample_counts[sid][type]:
                    all_keys[k] = True

            sorted_keys = sorted(all_keys.keys())
                    
            for k in sorted_keys:
                row = [k]
                for sid in sample_ids:
                    sc = sample_counts[sid][type]
                    ct = sc[k] if k in sc else 0
                    row.append(str(ct))
                cfh.write("\t".join(row) + "\n")

            cfh.write("\n\n")

        # ------------------------------------------------------
        # combine samples
        # ------------------------------------------------------
        find_unique_MEs(cfh, sample_inds, args.output_dir)
    
            
if __name__ == '__main__':
    main()


