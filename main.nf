#!/usr/bin/env nextflow
/*
 * Defines the pipeline input parameters (with a default value for each one).
 * Each of the following parameters can be specified as command line options.
 */
 params.transcript_gtf = "/data1/shahs3/users/preskaa/APS010.1_Archive/bambu_out/multiSample_NDR_0.1/detected_transcripts.gtf"
 params.transdecoder_gff3 = "/data1/shahs3/users/preskaa/APS010.1_Archive/transdecoder/NDR_0.1/transcripts.fa.transdecoder.gff3"
 params.transcripts_fasta = "/data1/shahs3/users/preskaa/APS010.1_Archive/gffread/NDR_0.1/transcripts.fa"
 params.outdir = "/data1/shahs3/users/preskaa/A673/data/APS010.1_PG3_v_Swissprot/transdecoder_ORF_viz" 
 params.sample_id = "A673"

 workflow {
    GTF2BED(params.transcript_gtf)
    GTF2GFF3(params.transcript_gtf)
    GENOME_ALIGNED_GFF3(GTF2GFF3.out, params.transdecoder_gff3, params.transcripts_fasta)
    GFF3_TO_BED(GENOME_ALIGNED_GFF3.out)
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
   publishDir "${params.outdir}/cds_bed/${params.sample_id}"
   container "quay.io/preskaa/biopython:v250221"

   input:
   path "transcripts.bed"

   output:
   path "transcripts.fasta.transdecoder.genome.bed"

   script:
   """
   #!/usr/bin/env python
   import sys
   import pandas as pd
   tx_bed = pd.read_csv("transcripts.bed", sep="\t", skiprows=1, names=col_names)
   tx_bed = tx_bed.drop_duplicates()
   tx_bed.to_csv("transcripts.fasta.transdecoder.genome.bed", 
                  sep="\t", 
                  header=False, 
                  index=False)
   # edit tag line
   track_line = 'track name="detected transcripts" visibility=2 itemRgb="On"\n'

   with open("transcripts.fasta.transdecoder.genome.bed", "r") as original:
      data = original.read()

   with open("transcripts.fasta.transdecoder.genome.bed", "w") as modified:
      modified.write(track_line + data)
   """

}


