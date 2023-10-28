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
def filter_and_index_csv_file(csv_dir, output_dir, output_suffix, cfile, filters):
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

                # -------------------------------------------------
                # write MEI to output file and add to index
                # -------------------------------------------------
                if not filter:
                    ofh.write(",".join(row) + "\n")
                    n_written += 1
                    # index
                    # TODO - could consider adding length and/or insertion sequence depending on our definition of "the same"
                    key = ":".join([chrom, pos, strand, iseq])
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
def find_unique_MEs(ME_inds):
    sample_ids = sorted(ME_inds.keys())
    
    # map each ME to the samples in which it appears
    unique_MEs = {}
    
    for sid in sample_ids:
        s_ind = ME_inds[sid]
        for key in s_ind:
            if key not in unique_MEs:
                unique_MEs[key] = []
            unique_MEs[key].append(sid)

    # sample count histogram
    sc_hist = {}
    for k in unique_MEs:
        sc = len(unique_MEs[k])
        if sc not in sc_hist:
            sc_hist[sc] = 0
        sc_hist[sc] += 1
            
    print("num unique ME(s) = " + str(len(unique_MEs)))
    keys = [k for k in sc_hist.keys()]
    for k in sorted(keys, key=lambda x: int(x), reverse=True):
        print(str(k) + " : " + str(sc_hist[k]))
    
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
    # SVA filters
    parser.add_argument('--sva_excluded_repeat_types', required=False, help='Exclude/filter SVAs whose overlapping repeat type is in this list.')
    # Alu filters
    parser.add_argument('--alu_excluded_repeat_types', required=False, help='Exclude/filter Alus whose overlapping repeat type is in this list.')
    # LINE filters
    parser.add_argument('--line_excluded_repeat_types', required=False, help='Exclude/filter LINEs whose overlapping repeat type is in this list.')
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

    filters = {
        'SVA': { 'repeat_types': ex_sva_rep_types_d },
        'ALU' : { 'repeat_types': ex_alu_rep_types_d },
        'LINE1': { 'repeat_types': ex_line_rep_types_d }
    }
        
    # read, filter, and index csv_files
    csv_files = read_csv_dir(args.csv_dir)
    sample_inds = {}
    for cf in csv_files:
        (sample_id, index) = filter_and_index_csv_file(args.csv_dir, args.output_dir, args.output_suffix, cf, filters)
        sample_inds[sample_id] = index

    # ------------------------------------------------------
    # write summary/counts file
    # ------------------------------------------------------
    counts_file = "summary-counts-" + args.output_suffix + ".tsv"
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
    find_unique_MEs(sample_inds)
    
            
if __name__ == '__main__':
    main()


