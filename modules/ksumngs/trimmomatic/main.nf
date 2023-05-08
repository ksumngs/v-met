process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::trimmomatic=0.39" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2':
        'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readtype = meta.single_end ? 'SE' : 'PE'
    def trimmed = meta.single_end ?
        "${prefix}_trimmed.fastq.gz" :
        "${prefix}_trimmed_R1.fastq.gz /dev/null ${prefix}_trimmed_R2.fastq.gz /dev/null"
    def jmemstring = task.memory.toMega() + 'M'
    """
    trimmomatic \\
            -Xmx${jmemstring} \\
            ${readtype} \\
            -threads ${task.cpus} \\
            ${reads} \\
            ${trimmed} \\
            ${args} \\
        2> >(tee ${prefix}.log >&2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic -version)
    END_VERSIONS
    """
}
