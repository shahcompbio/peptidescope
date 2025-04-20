// create igv-report using proteins of interest
process IGVREPORT {
    tag 'process_low'
    publishDir "${params.outdir}/visualization/${params.sample_id}", mode: 'copy', overwrite: true

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/igv-reports:1.14.1--pyh7e72e81_0'
        : 'biocontainers/igv-reports:1.14.1--pyh7e72e81_0'}"

    input:
    path regions_bed
    // regions of interest
    path refgenome
    // ref genome sample was aligned to
    tuple path(bigwig), path(peptides), path(transcripts)

    output:
    path "peptidescope_report.html", emit: html

    script:
    """
   create_report ${regions_bed} \
      --fasta ${refgenome} \
      --genome hg38 \
      --info-columns ID \
      --zero_based true \
      --tracks ${bigwig} ${peptides} ${transcripts} \
      --output peptidescope_report.html
   """
}
