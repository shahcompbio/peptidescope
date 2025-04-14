#!/usr/bin/env nextflow
/*
 * Defines the pipeline input parameters (with a default value for each one).
 * Each of the following parameters can be specified as command line options.
 */
 params.transcript_gtf = "/data1/shahs3/users/preskaa/APS010.1_Archive/bambu_out/multiSample_NDR_0.1/detected_transcripts.gtf"
 params.transdecoder_gff3 = "/data1/shahs3/users/preskaa/APS010.1_Archive/transdecoder/NDR_0.1/transcripts.fa.transdecoder.gff3"
 params.transcripts_fasta = "/data1/shahs3/users/preskaa/APS010.1_Archive/gffread/NDR_0.1/transcripts.fa"
 params.out_dir = "/data1/shahs3/users/preskaa/A673/data/APS010.1_PG3_v_Swissprot/transdecoder_ORF_viz" 
 params.sample = "A673"

 workflow {
    GTF2BED(params.transcript_gtf)
    GTF2GFF3(params.transcript_gtf)
 }
// convert to gtf to bed file
 process GTF2BED {
    tag "GTF2BED"
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
    tag "GTF2GFF3"
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
    tag "GENOMEGFF3"
    container "quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"
    input:
    path 'transcripts.gff3'
    path transdecoder_gff3
    path transcripts_fasta

    output:
    path "transcripts.transdecoder.genome.gff3"

    """
    perl scripts/cdna_alignment_orf_to_genome_orf.pl
    """
 }
