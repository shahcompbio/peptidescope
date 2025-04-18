#!/usr/bin/env nextflow
/*
 * Defines the pipeline input parameters (with a default value for each one).
 * Each of the following parameters can be specified as command line options.
 */
 params.transcript_gtf = "/data1/shahs3/users/preskaa/APS010.1_Archive/bambu_out/multiSample_NDR_0.1/detected_transcripts.gtf"
 params.transdecoder_gff3 = "/data1/shahs3/users/preskaa/APS010.1_Archive/transdecoder/NDR_0.1/transcripts.fa.transdecoder.gff3"
 params.transcripts_fasta = "/data1/shahs3/users/preskaa/APS010.1_Archive/gffread/NDR_0.1/transcripts.fa"
 params.outdir = "/data1/shahs3/users/preskaa/A673/data/APS010.1_PG3_v_Swissprot/transdecoder_ORF_viz"
 params.fragpipe_dir = "/data1/shahs3/users/preskaa/APS010.1_Archive/fragpipe/spike_in"
 params.sample_id = "A673"
 params.igv_report = false
 params.bigwig = "/data1/shahs3/isabl_data_lake/analyses/24/35/42435/results/minimap2/bigwig/SHAH_H003599_T01_01_TR01_R1.bedGraph"
 params.protein_list = "/data1/shahs3/users/preskaa/A673/data/APS010.1_PG3_v_Swissprot/peptidescope/noblastmatch_proteins.txt"
 params.ref_genome = "/data1/shahs3/reference/ref-sarcoma/GRCh38/v45/GRCh38.primary_assembly.genome.fa"
 params.ref_genome_index = "/data1/shahs3/reference/ref-sarcoma/GRCh38/v45/GRCh38.primary_assembly.genome.fa.fai"

 // modules
 include { BEDTOOLS_SLOP } from 'modules/nf-core/bedtools/slop/main'
 include { BEDTOOLS_INTERSECT } from 'modules/nf-core/bedtools/intersect/main'
 workflow {
    GTF2BED(params.transcript_gtf)
    GTF2GFF3(params.transcript_gtf)
    GENOME_ALIGNED_GFF3(GTF2GFF3.out, params.transdecoder_gff3, params.transcripts_fasta)
    GFF3_TO_BED(GENOME_ALIGNED_GFF3.out)
    transcripts_bed = DUPS(GFF3_TO_BED.out)
    peptides_bed = PEPTIDE_BED(params.fragpipe_dir, transcripts_bed)
    if (params.igv_report) {
      meta = [id: params.sample_id, description: "igv_report"]
      regions_bed = INTERESTING_BED(params.protein_list, transcripts_bed)
      GENOME_SIZES(params.ref_genome_index)
      expanded_regions_bed = BEDTOOLS_SLOP(meta, regions_bed, GENOME_SIZES.out){ext.args = "-b 1000"}
      // intersect tracks of interest
      filtered_transcripts = BEDTOOLS_INTERSECT(meta, transcripts_bed, expanded_regions_bed, GENOME_SIZES.out)
      filtered_peptides = BEDTOOLS_INTERSECT(meta, peptides_bed, expanded_regions_bed, GENOME_SIZES.out)
      filtered_bigwig = BEDTOOLS_INTERSECT(meta, params.bigwig, expanded_regions_bed, GENOME_SIZES.out)
      //create igv report
      IGV_REPORT(regions_bed, params.ref_genome, filtered_bigwig, filtered_peptides, filtered_transcripts)

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
process GFF3_TO_BED{
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
process DUPS{
   tag "process_low"
   publishDir "${params.outdir}/visualization/${params.sample_id}"
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
    publishDir "${params.outdir}/visualization/${params.sample_id}"

    input:
    path fragpipe_dir
    path transcripts_genome_bed // transcriptome bed file in ref genome coordinates

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
process INTERESTING_BED{
   tag "process_low"
   container "quay.io/preskaa/biopython:v250221"

   input:
   path "proteins.txt" // list of proteins of interest
   path "transcripts.genome.bed" // transcript bed file

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
process GENOME_SIZES{
   tag "process_low"
   container "quay.io/nf-core/ubuntu:22.04"

   input:
   path genome_index // ref genome index

   output:
   path "genome_sizes.txt"

   script:
   """
   cut -f1,2 ${genome_index} > genome_sizes.txt
   """
}

// generate igv report
process IGV_REPORT{
   tag "process_low"
   container "quay.io/biocontainers/igv-reports:1.14.1--pyh7e72e81_0"
   publishDir "${params.outdir}/visualization/${params.sample_id}"

   input:
   path regions_bed // regions of interest
   path refgenome // ref genome sample was aligned to
   path filtered_bigwig // filtered bigwig
   path filtered_peptides_bed // peptides filtered
   path filtered_transcripts_bed // filtered transcripts bed file

   output:
   path "peptides_igv_report.html"

   script:
   """
   create_report ${regions_bed} \
      --fasta ${refgenome} \
      --genome hg38 \
      --info-columns ID \
      --zero_based true \
      --tracks ${filtered_bigwig} ${filtered_peptides_bed} ${filtered_transcripts_bed} \
      --output peptides_igv_report.html
   """
}
