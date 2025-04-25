#!/usr/bin/env nextflow


// modules
include { BEDTOOLS_SLOP      } from './modules/local/bedtools/slop/main'
include { BEDTOOLS_INTERSECT } from './modules/local/bedtools/intersect/main'
include { IGVREPORT          } from './modules/local/igvreport/main'

workflow {
   GTF2BED(params.transcript_gtf)
   GTF2GFF3(params.transcript_gtf)
   GENOME_ALIGNED_GFF3(GTF2GFF3.out, params.transdecoder_gff3, params.transcripts_fasta)
   GFF3_TO_BED(GENOME_ALIGNED_GFF3.out)
   transcripts_bed = DUPS(GFF3_TO_BED.out)
   peptides_bed = PEPTIDE_BED(params.fragpipe_dir, transcripts_bed)
   if (params.igv_report) {
      regions_bed = INTERESTING_BED(params.protein_list, transcripts_bed)
      GENOME_SIZES(params.ref_genome_index)
      expanded_regions_bed = BEDTOOLS_SLOP(regions_bed, GENOME_SIZES.out)
      // intersect tracks with region of interest; here i pay the price of lack
      // of metamaps ...
      bw_ch = channel.of([1, file(params.bigwig)])
      peptides_ch = peptides_bed.map { f -> [2, file(f)] }
      transcripts_ch = transcripts_bed.map { f -> [3, file(f)] }
      tracks_ch = bw_ch.concat(peptides_ch, transcripts_ch)
      // intersect tracks of interest
      combined_ch = tracks_ch.combine(expanded_regions_bed)

      igv_tracks = BEDTOOLS_INTERSECT(combined_ch).toSortedList { a, b -> a[0] <=> b[0] }.map { v -> v.collect { it[1] } }

      IGVREPORT(regions_bed, params.ref_genome, igv_tracks)
   }
}
// convert to gtf to bed file
process GTF2BED {
   tag "process_low"
   container "quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"

   input:
   path transcript_gtf

   output:
   path "transcripts.bed"

   script:
   """
    gtf_to_bed.pl ${transcript_gtf} > transcripts.bed
    """
}
// convert the predictions to a genome based bed file as well
process GTF2GFF3 {
   tag "process_low"
   container "quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"

   input:
   path transcript_gtf

   output:
   path "transcripts.gff3"

   script:
   """
    gtf_to_alignment_gff3.pl ${transcript_gtf} > transcripts.gff3
    """
}
// generate genome coordinate based CDS annotation
process GENOME_ALIGNED_GFF3 {
   tag "process_low"
   container "quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"

   input:
   path transcripts_gff3
   path transdecoder_gff3
   path transcripts_fasta

   output:
   path "transcripts.transdecoder.genome.gff3"

   script:
   """
    transcript_orf_to_genome_orf.pl ${transdecoder_gff3} ${transcripts_gff3} \
    ${transcripts_fasta} > transcripts.transdecoder.genome.gff3
    """
}
// generate bed file
process GFF3_TO_BED {
   tag "process_low"
   container "quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"

   input:
   path genome_gff3

   output:
   path "transcripts.fasta.transdecoder.genome.bed"

   script:
   """
   gff3_file_to_bed.pl ${genome_gff3} > transcripts.fasta.transdecoder.genome.bed
   """
}
// remove duplicate entries (this is because we stink at perl)
process DUPS {
   tag "process_low"
   publishDir "${params.outdir}/visualization/${params.sample_id}", mode: 'copy', overwrite: true
   container "quay.io/preskaa/biopython:v250221"

   input:
   path "transcripts.bed"

   output:
   path "transcripts.fasta.transdecoder.genome.bed"

   script:
   """
#!/usr/bin/env python
import pandas as pd
col_names = ["chrom", "chromStart", "chromEnd", "name", 
               "score", "strand", "thickStart", "thickEnd",
               "itemRgb", "blockCount", "blockSizes", "blockStarts"]
tx_bed = pd.read_csv("transcripts.bed", sep="\t", skiprows=1, names=col_names)
tx_bed = tx_bed.drop_duplicates()
tx_bed.to_csv("transcripts.fasta.transdecoder.genome.bed", 
               sep="\t", 
               header=False, 
               index=False)
# edit tag line
track_line = 'track name="detected transcripts" visibility=2 itemRgb="On"\\n'

with open("transcripts.fasta.transdecoder.genome.bed", "r") as original:
   data = original.read()

with open("transcripts.fasta.transdecoder.genome.bed", "w") as modified:
   modified.write(track_line + data)
"""
}

// generate bed file with peptides
process PEPTIDE_BED {
   tag "process_low"
   container "quay.io/preskaa/biopython:v250221"
   publishDir "${params.outdir}/visualization/${params.sample_id}", mode: 'copy', overwrite: true

   input:
   path fragpipe_dir
   path transcripts_genome_bed

   output:
   path "peptides.bed"

   script:
   """
    peptide_bedfile.py ${fragpipe_dir} ${transcripts_genome_bed} peptides.bed
    """
}

/*
 * will refactor the following into a subworkflow called IGV_REPORT
 */

// generate bedfile of proteins of interest 
process INTERESTING_BED {
   tag "process_low"
   container "quay.io/preskaa/biopython:v250221"
   publishDir "${params.outdir}/visualization/${params.sample_id}", mode: 'copy', overwrite: true

   input:
   path "proteins.txt"
   // list of proteins of interest
   path "transcripts.genome.bed"

   output:
   path "igv_report.proteins.bed"

   script:
   """
#!/usr/bin/env python
import pandas as pd
# read in protein list
with open("proteins.txt", 'r') as file:
    proteins = [line.strip() for line in file]
col_names = ["chrom", "chromStart", "chromEnd", "name", 
               "score", "strand", "thickStart", "thickEnd",
               "itemRgb", "blockCount", "blockSizes", "blockStarts"]
tx_bed = pd.read_csv("transcripts.genome.bed", sep="\t", skiprows=1, names=col_names)
tx_bed = tx_bed[tx_bed["name"].apply(
    lambda x: any(f"{prot}.p" in x for prot in proteins)
)]
tx_bed.to_csv("igv_report.proteins.bed", 
               sep="\t", 
               header=False, 
               index=False)
"""
}

// generate genome sizes for bedtools slop
process GENOME_SIZES {
   tag "process_low"
   container "quay.io/nf-core/ubuntu:22.04"

   input:
   path genome_index

   output:
   path "genome_sizes.txt"

   script:
   """
   cut -f1,2 ${genome_index} > genome_sizes.txt
   """
}
