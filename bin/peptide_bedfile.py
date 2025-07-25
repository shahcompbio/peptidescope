#!/usr/bin/env python
import numpy as np
import pandas as pd
import os
import sys

#paths
#TODO: there is an occasional bug where the peptides appear to be one AA short; need to have code review
archive = sys.argv[1]
tx_bed_path = sys.argv[2]
pep_bed_path = sys.argv[3]

# useful functions
def get_subdirectories(folder_path):
  """
  Returns a list of all subdirectories within the given folder path.

  Args:
    folder_path: The path to the folder to search.

  Returns:
    A list of strings, where each string is the full path to a subdirectory.
    Returns an empty list if the folder does not exist or has no subdirectories.
  """
  if not os.path.isdir(folder_path):
    return []

  subdirectories = [f.path for f in os.scandir(folder_path) if f.is_dir()]
  return subdirectories

def pep_ref_pos(ORF_start_idx, ORF_start, block_ends, 
                pep_pos, block_sizes, block_starts, ORF_end):
    """
    determine peptide genomic coordinates
    ORF_start_idx: index of exon in which ORF begins
    ORF_start: genomic position of ORF start
    block_ends: genomic coordinates of exon ends
    pep_pos: cDNA position in protein coord space where peptide starts, accounting for strand
    block_sizes: sizes of exons
    """
    # reset interval
    interval = ORF_start_idx
    # determine peptide start position in genomic coordinates
    block_size = block_ends[interval]-ORF_start
    bps = pep_pos
    while bps >= block_size and interval + 1 < len(block_sizes):
        # subtract off cDNA from block
        bps = bps-block_size
        # jump to next block
        interval = interval + 1
        # fetch block size of next block
        block_size = block_sizes[interval]
    # set genomic peptide start position
    if ORF_start_idx == interval:
        genomic_pep_pos = bps+ORF_start
    # fresh outta peptides
    elif bps == 0:
        genomic_pep_pos = ORF_end
    else:
        genomic_pep_pos = bps+block_starts[interval]
    assert genomic_pep_pos <= ORF_end, \
        f"Peptide position {genomic_pep_pos} exceeds ORF end {ORF_end}"
    return genomic_pep_pos

def calc_orf_size(block_sizes, start_interval, end_interval, 
                  block_starts, block_ends, ORF_start, ORF_end):
    """
    compute ORF cDNA length
    """
    # check if ORF is in single block; will also catch single exons
    if start_interval == end_interval:
        ORF_size = ORF_end-ORF_start
    else:
        ORF_size = 0
        for i in np.arange(start_interval, end_interval + 1):
            if i == start_interval:
                block_size = block_ends[i] - ORF_start
            elif i == end_interval:
                block_size = ORF_end - block_starts[i]
            else:
                block_size = block_sizes[i]
            ORF_size += block_size
    return ORF_size

# get subdirectories for enzymes
enzyme_subdirs = get_subdirectories(archive)
detected_df1 = pd.DataFrame()
for i in np.arange(0, len(enzyme_subdirs)):
    enzyme_dir = enzyme_subdirs[i]
    fa_path = os.path.join(enzyme_dir, "peptide.tsv")
    if not os.path.exists(fa_path):
        print(f"peptide.tsv not found in {enzyme_dir}")
        continue
    temp = pd.read_csv(fa_path, sep="\t")
    # protein csv
    fa_path = os.path.join(enzyme_dir, "protein.tsv")
    protein_df = pd.read_csv(fa_path, sep="\t")
    protein_df = protein_df[["Protein", "Length"]]
    protein_df.columns = ["Protein", "Protein Length"]
    temp = pd.merge(temp, protein_df, on="Protein", how="left")
    detected_df1 = pd.concat([detected_df1, temp])


# drop indistinguishable proteoforms
unique_df = detected_df1[detected_df1["Mapped Proteins"].isna()]
# load in transcript structures
col_names = ["chrom", "chromStart", "chromEnd", "name", 
                  "score", "strand", "thickStart", "thickEnd",
                  "itemRgb", "blockCount", "blockSizes", "blockStarts"]
tx_bed = pd.read_csv(tx_bed_path, sep="\t", skiprows=1, names=col_names)
# just in case but should be resolved now
tx_bed = tx_bed.drop_duplicates()
## write our peptide bedfile
# group on protein ids
# drop duplicates for peptides which detected by multiple enzymes
unique_df = unique_df.drop_duplicates(subset=["Peptide", "Protein Start", "Protein End", "Protein"])
protein_groups = unique_df.groupby(by="Protein")
# initiate dataframe
data = []
# loop through protein groups; could be faster but should work
for protein, group in protein_groups:
    # fetch transcript info
    tx_id = protein
    bedrow = tx_bed[tx_bed["name"].str.contains(f"{tx_id}.p")]
    if len(bedrow) < 1:
        print(f"no transcript match for {tx_id}")
        continue
    assert len(bedrow) == 1, f"Expected 1 match for {tx_id}, found {len(bedrow)}"
    bedrow = bedrow.squeeze()
    block_sizes = [int(x) for x in bedrow["blockSizes"].split(",")]
    block_starts = [int(x) for x in bedrow["blockStarts"].split(",")]
    ORF_start = bedrow["thickStart"]
    ORF_end = bedrow["thickEnd"]
    exon_count = bedrow["blockCount"]
    # get block starts in chromosomal coordinates
    block_starts = bedrow["chromStart"] + np.array(block_starts)
    # get block ends in chromosomal coordinates
    block_ends = block_starts+np.array(block_sizes)
    # determine which block the ORF start is in
    interval_idx = np.where((block_starts <= ORF_start) & (ORF_start <= block_ends))[0]
    start_interval = interval_idx[0]
    # determine block in which ORF ends ...
    end_idx = np.where((block_starts <= ORF_end) & (ORF_end <= block_ends))[0]
    end_interval = end_idx[0]
    # now get size of ORF in cDNA
    ORF_size = calc_orf_size(block_sizes, start_interval, end_interval,
                             block_starts, block_ends, ORF_start, ORF_end)
    i = 1
    # determine peptide positions
    for _, row in group.iterrows():
        peplen = 3*(row["Peptide Length"])
        if bedrow["strand"] == "+":
            # get peptide start and end; correcting for differences in indexing
            pepstart = 3*(row["Protein Start"]-1)
            pepend = 3*(row["Protein End"]-1)
            # determine peptide start position in genomic coordinates
            genomic_pepstart = pep_ref_pos(start_interval, ORF_start, block_ends, 
                                           pepstart, block_sizes, block_starts, ORF_end)
            # determine peptide end position in genomic coordinates
            genomic_pepend = pep_ref_pos(start_interval, ORF_start, block_ends, 
                                         pepend, block_sizes, block_starts, ORF_end)
        else:
            # negative strand
            # flip protein start end coordinates
            pepstart = ORF_size - 3*(row["Protein End"])
            pepend = ORF_size - 3*(row["Protein Start"]-1)
            # determine peptide start position in genomic coordinates
            genomic_pepstart = pep_ref_pos(start_interval, ORF_start, block_ends, 
                                           pepstart, block_sizes, block_starts, ORF_end)
            # determine peptide end position in genomic coordinates
            genomic_pepend = pep_ref_pos(start_interval, ORF_start, block_ends, 
                                         pepend, block_sizes, block_starts, ORF_end) 
        # check that coordinates make sense
        assert ORF_end >= genomic_pepend, \
            f"peptide ends after ORF {genomic_pepend} > {ORF_end} for {tx_id}"
        assert ORF_start <= genomic_pepstart, \
            f"peptide starts before ORF {genomic_pepstart} < {ORF_start} for {tx_id}"
        data.append({
            'chrom':bedrow['chrom'],
            'chromStart':bedrow['chromStart'],
            'chromEnd':bedrow['chromEnd'],
            'name': f"{tx_id}_peptide_{i}",
            'score': bedrow['score'],
            'strand': bedrow['strand'],
            'thickStart':genomic_pepstart,
            'thickEnd': genomic_pepend,
            'itemRgb': '0',
            'blockCount': bedrow["blockCount"],
            'blockSizes': bedrow['blockSizes'],
            'blockStarts': bedrow['blockStarts']
        })
        i += 1
test_pep_bed = pd.DataFrame(data)
# make our test bed file
test_pep_bed.to_csv(pep_bed_path, sep="\t", header=False, index=False)
track_line = 'track name="unique peptides" description="detected peptides" visibility=2 itemRgb="On"\n'

with open(pep_bed_path, "r") as original:
    data = original.read()

with open(pep_bed_path, "w") as modified:
    modified.write(track_line + data)