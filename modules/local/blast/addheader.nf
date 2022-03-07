process BLAST_ADDHEADER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.1.0_cv2/biocontainers_v1.1.0_cv2.img' :
        'docker.io/biocontainers/biocontainers:v1.1.0_cv2' }"

    input:
    tuple val(meta), path(blast, stageAs: 'blast.txt')
    path(header, stageAs: 'header.txt')

    output:
    tuple val(meta), path("*.blast.tsv"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${header} ${blast} > ${prefix}.blast.tsv
    """
}
