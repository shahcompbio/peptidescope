// filter input files to only those overlapping the expanded regions

process BEDTOOLS_INTERSECT {
    tag 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--h13024bc_3':
        'biocontainers/bedtools:2.31.1--h13024bc_3' }"

    input:
    tuple val(id), path(track), path(region)

    output:
    tuple val(id), path("filtered.*")

    script:
    def track_name = track.getName()
    """
    bedtools intersect -a ${track} -b ${region} > filtered.${track_name}
    """
}
