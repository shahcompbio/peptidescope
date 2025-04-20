// expand regions by 1000 bp on either side

process BEDTOOLS_SLOP {
    tag 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--h13024bc_3':
        'biocontainers/bedtools:2.31.1--h13024bc_3' }"

    input:
    path raw_regions_bed
    path genome_sizes

    output:
    path "expanded.regions.bed", emit: bed

    script:
    """
    bedtools slop -i ${raw_regions_bed} -g ${genome_sizes} -b 100000 > expanded.regions.bed
    """
}
