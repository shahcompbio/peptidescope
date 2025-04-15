 #!/usr/bin/env python
from tqdm import tqdm
import numpy as np
import pandas as pd
import os
from gtfparse import read_gtf
import sys

#paths

archive = sys.argv[1]
gtfpath = "/data1/shahs3/users/preskaa/APS010.1_Archive/bambu_out/multiSample_NDR_0.1/detected_transcripts.gtf"
tx_bed_path = "/data1/shahs3/users/preskaa/A673/data/APS010.1_PG3_v_Swissprot/transdecoder_ORF_viz/visualization/A673/transcripts.fasta.transdecoder.genome.bed"
test_pep_bed_path = "/data1/shahs3/users/preskaa/A673/data/APS010.1_PG3_v_Swissprot/transdecoder_ORF_viz/visualization/A673/test_peptides.bed"

# useful functions
def pep_ref_pos(ORF_start_idx, ORF_start, block_ends, 
                pep_pos, block_sizes, block_starts):
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
    while bps > block_size:
        # subtract off cDNA from block
        bps = bps-block_size
        # jump to next block
        interval = interval + 1
        # fetch block size of next block
        block_size = block_sizes[interval]
    # set genomic peptide start position
    if ORF_start_idx == interval:
        genomic_pepstart = bps+ORF_start
    else:
        genomic_pepstart = bps+block_starts[interval]
    return genomic_pepstart

# load enzymes (will refactor later for general use)
enzymes = ["argc", "aspn", "gluc", 
           "in-house_chymotrypsin", "lysc", "lysn", 
           "proalanase", "trypsin"]
detected_df1 = pd.DataFrame()
for i in tqdm(np.arange(0, len(enzymes))):
    enzyme = enzymes[i]
    fa_path = os.path.join(archive, enzyme + "_diaPASEF_groupFDR", "peptide.tsv")
    temp = pd.read_csv(fa_path, sep="\t")
    # protein csv
    fa_path = os.path.join(archive, enzyme + "_diaPASEF_groupFDR", "protein.tsv")
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
protein_groups = unique_df.groupby(by="Protein")
# initiate dataframe
data = []
# loop through protein groups; could be faster but should work
for protein, group in protein_groups:
    # fetch transcript info
    tx_id = protein
    bedrow = tx_bed[tx_bed["name"].str.contains(tx_id)]
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
    interval_idx = np.where((block_starts <= ORF_start) & (ORF_start < block_ends))[0]
    start_interval = interval_idx[0]
    # determine block in which ORF ends ...
    end_idx = np.where((block_starts <= ORF_end) & (ORF_end < block_ends))[0]
    end_interval = end_idx[0]
    # now get size of ORF in cDNA
    ORF_size = 0
    for i in np.arange(0, len(block_sizes)):
        if i == start_interval:
            block_size = block_ends[i] - ORF_start
        elif i == end_interval:
            block_size = ORF_end - block_starts[i]
        else:
            block_size = block_sizes[i]
        ORF_size += block_size
    i = 1
    # determine peptide positions
    for _, row in group.iterrows():
        peplen = 3*(row["Peptide Length"])
        if bedrow["strand"] == "+":
            # get peptide start and end
            pepstart = 3*(row["Protein Start"]-1)
            pepend = 3*(row["Protein End"])
            # determine peptide start position in genomic coordinates
            genomic_pepstart = pep_ref_pos(start_interval, ORF_start, block_ends, 
                                           pepstart, block_sizes, block_starts)
            # determine peptide end position in genomic coordinates
            genomic_pepend = pep_ref_pos(start_interval, ORF_start, block_ends, 
                                         pepend, block_sizes, block_starts)
        else:
            # negative strand
            # flip protein start end coordinates
            pepstart = ORF_size - 3*(row["Protein End"])
            pepend = ORF_size - 3*(row["Protein Start"]-1)
            # determine peptide start position in genomic coordinates
            genomic_pepstart = pep_ref_pos(start_interval, ORF_start, block_ends, 
                                           pepstart, block_sizes, block_starts)
            # determine peptide end position in genomic coordinates
            genomic_pepend = pep_ref_pos(start_interval, ORF_start, block_ends, 
                                         pepend, block_sizes, block_starts)
            # break
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
            'blockStarts': bedrow['blockStarts'],
            
        })
        i += 1
test_pep_bed = pd.DataFrame(data)
# make our test bed file
test_pep_bed.to_csv(test_pep_bed_path, sep="\t", header=False, index=False)
track_line = 'track name="unique peptides" description="detected peptides" visibility=2 itemRgb="On"\n'

with open(test_pep_bed_path, "r") as original:
    data = original.read()

with open(test_pep_bed_path, "w") as modified:
    modified.write(track_line + data)