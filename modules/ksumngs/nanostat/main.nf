process NANOSTAT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::nanostat=1.6.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanostat:1.6.0--pyhdfd78af_0' :
        'quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_NanoStats"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Figure out what kind of file we are collecting stats on
    // Start with assuming this is a summary file, as the extension
    // is less consistent
    def analysis_flag = 'summary'
    analysis_flag = (reads.getName().contains('.fa') || reads.getName().contains('.fasta')) ? 'fasta' : analysis_flag
    analysis_flag = (reads.getName().contains('.fq') || reads.getName().contains('.fastq')) ? 'fastq' : analysis_flag
    analysis_flag = (reads.getName().contains('.bam')) ? 'bam' : analysis_flag
    """
    NanoStat \\
            -t ${task.cpus} \\
            --${analysis_flag} ${reads} \\
            ${args} \\
        > ${prefix}_NanoStats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanostat: \$(NanoStat -v | sed 's/NanoStat //')
    END_VERSIONS
    """
}
