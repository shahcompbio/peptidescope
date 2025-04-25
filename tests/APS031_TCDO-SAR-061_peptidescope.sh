#!/bin/bash
## activate nf-core conda environment
source /home/preskaa/miniforge3/bin/activate nf-core

module load java/20.0.1
## specify params
outdir=/data1/shahs3/users/preskaa/ThreeByThreeSarcoma/data/APS031_3x3_proteomics/peptidescope
pipelinedir=$HOME/peptidescope
transcript_gtf=/data1/shahs3/users/preskaa/APS022.1_Archive/ont_rna/bambu_out/multiSample_NDR_0.1/detected_transcripts.gtf
transcripts_fasta=/data1/shahs3/users/preskaa/APS022.1_Archive/ont_rna/gffread/NDR_0.1/transcripts.fa
transdecoder_gff3=/data1/shahs3/users/preskaa/APS022.1_Archive/ont_rna/transdecoder/NDR_0.1/transcripts.fa.transdecoder.gff3
sample_id=TCDO-SAR-061
fragpipe_dir=/data1/shahs3/users/preskaa/ThreeByThreeSarcoma/data/APS031_3x3_proteomics/fragpipe_out/TCDO_SAR_061/ONT
protein_list=${outdir}/${sample_id}_proteins.txt
bigwig=/data1/shahs3/isabl_data_lake/analyses/33/30/43330/results/minimap2/bigwig/SHAH_H003842_T01_01_TR01_R1.bedGraph
ref_genome=/data1/shahs3/reference/ref-sarcoma/GRCh38/v45/GRCh38.primary_assembly.genome.fa

mkdir -p ${outdir}
cd ${outdir}

nextflow run ${pipelinedir}/main.nf \
    -profile singularity,slurm \
    -work-dir ${outdir}/work \
    --outdir ${outdir} \
    --transcript_gtf ${transcript_gtf} \
    --transdecoder_gff3 ${transdecoder_gff3} \
    --transcripts_fasta ${transcripts_fasta} \
    --sample_id ${sample_id} \
    --fragpipe_dir ${fragpipe_dir} \
    --igv_report \
    --bigwig ${bigwig} \
    --protein_list ${protein_list} \
    --ref_genome ${ref_genome} \
    --ref_genome_index ${ref_genome}.fai \
    -resume