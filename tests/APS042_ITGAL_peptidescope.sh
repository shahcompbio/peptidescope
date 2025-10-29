#!/bin/bash
outdir=/data1/shahs3/users/preskaa/AMLproteogenomics/data/APS042_ITGAL_protein/peptidescope_out
transcript_gtf=/data1/shahs3/isabl_data_lake/analyses/53/26/45326/bambu_out/transcriptome_NDR_0.1/detected_transcripts.gtf
transdecoder_gff3=/data1/shahs3/isabl_data_lake/analyses/53/41/45341/results/transdecoder/SHAH_H003942_T01_01_TR01/transcripts.fa.transdecoder.gff3
transcripts_fasta=/data1/shahs3/isabl_data_lake/analyses/53/41/45341/results/gffread/SHAH_H003942_T01_01_TR01/transcripts.fa
sample_id=MV4-11
fragpipe_dir=/data1/shahs3/users/preskaa/AMLproteogenomics/data/APS042_ITGAL_protein/MV4-11
bigwig=/data1/shahs3/isabl_data_lake/analyses/49/26/44926/results/minimap2/bigwig/SHAH_H003942_T01_01_TR01_R1.bedGraph
protein_list=/data1/shahs3/users/preskaa/AMLproteogenomics/data/APS042_ITGAL_protein/protein_list.txt
protein_info=/data1/shahs3/isabl_data_lake/analyses/53/41/45341/results/proteome/SHAH_H003942_T01_01_TR01/protein_transcript_info.tsv
ref_genome=/data1/shahs3/reference/ref-sarcoma/GRCh38/v45/GRCh38.primary_assembly.genome.fa
ref_genome_index=/data1/shahs3/reference/ref-sarcoma/GRCh38/v45/GRCh38.primary_assembly.genome.fa.fai
# run peptidescope on MV4-11 AML sample
nextflow run $HOME/peptidescope/main.nf \
    -profile singularity \
    --outdir ${outdir} \
    --transcript_gtf ${transcript_gtf} \
    --transdecoder_gff3 ${transdecoder_gff3} \
    --transcripts_fasta ${transcripts_fasta} \
    --sample_id ${sample_id} \
    --fragpipe_dir ${fragpipe_dir} \
    --bigwig ${bigwig} \
    --protein_list ${protein_list} \
    --protein_info ${protein_info} \
    --ref_genome ${ref_genome} \
    --ref_genome_index ${ref_genome_index} \
    --igv_report \
    -resume
