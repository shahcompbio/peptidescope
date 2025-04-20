#!/bin/bash

## activate nf-core conda environment
source /home/preskaa/miniforge3/bin/activate nf-core

module load java/20.0.1
## specify params
outdir=/data1/shahs3/users/preskaa/A673/data/APS010.1_PG3_v_Swissprot/peptidescope
pipelinedir=$HOME/peptidescope
transcript_gtf=/data1/shahs3/users/preskaa/APS010.1_Archive/bambu_out/multiSample_NDR_0.1/detected_transcripts.gtf
transdecoder_gff3=/data1/shahs3/users/preskaa/APS010.1_Archive/transdecoder/NDR_0.1/transcripts.fa.transdecoder.gff3
sample_id=A673
fragpipe_dir=/data1/shahs3/users/preskaa/APS010.1_Archive/fragpipe/spike_in
protein_list=${outdir}/noblast_match_proteins.txt
bigwig=/data1/shahs3/isabl_data_lake/analyses/24/35/42435/results/minimap2/bigwig/SHAH_H003599_T01_01_TR01_R1.bedGraph
ref_genome=/data1/shahs3/reference/ref-sarcoma/GRCh38/v45/GRCh38.primary_assembly.genome.fa

mkdir -p ${outdir}
cd ${outdir}

nextflow run ${pipelinedir}/main.nf \
    -profile singularity \
    -work-dir ${outdir}/work \
    --outdir ${outdir} \
    --transcript_gtf ${transcript_gtf} \
    --transdecoder_gff3 ${transdecoder_gff3} \
    --sample_id ${sample_id} \
    --fragpipe_dir ${fragpipe_dir} \
    --igv_report \
    --bigwig ${bigwig} \
    --protein_list ${protein_list} \
    --ref_genome ${ref_genome} \
    --ref_genome_idx ${ref_genome}.fai \
   # -resume
