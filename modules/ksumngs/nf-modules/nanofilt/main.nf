process NANOFILT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::nanofilt=2.8.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanofilt:2.8.0--py_0':
        'quay.io/biocontainers/nanofilt:2.8.0--py_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -cdf ${reads} \\
        | NanoFilt \\
            --logfile ${prefix}.nanofilt.log \\
            ${args} \\
        | gzip \\
        > ${prefix}_trimmed.fastq.gz

    # Workaround for no log file output
    # See https://github.com/wdecoster/nanofilt/issues/50
    echo " " >> ${prefix}.nanofilt.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanofilt: \$(NanoFilt -v | sed 's/NanoFilt //')
    END_VERSIONS
    """
}
