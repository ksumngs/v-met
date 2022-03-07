def VERSION = '1.2' // Version information not provided by tool on CLI

process KRAKENTOOLS_EXTRACT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::krakentools=1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakentools:1.2--pyh5e36f6f_0':
        'quay.io/biocontainers/krakentools:1.2--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(reads), path(krakenfile), path(kreport)
    val(taxids)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path 'versions.yml'                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readsFlag = (meta.single_end) ? "-s ${reads}" : "-s1 ${reads[0]} -s2 ${reads[1]}"
    def filteredFlag = (meta.single_end) ? "-o ${prefix}_filtered.fastq" : "-o ${prefix}_filtered_R1.fastq -o2 ${prefix}_filtered_R2.fastq"
    def kreportFlag = kreport ? "-r ${kreport}" : ''
    """
    # Kraken files should be gzipped, and are **VERY** large, so use a named pipe to
    # prevent storage overload
    mkfifo krakenfile
    gzip -cdf ${krakenfile} > krakenfile &
    extract_kraken_reads.py \\
        -k krakenfile \\
        ${readsFlag} \\
        ${filteredFlag} \\
        -t ${taxids} \\
        --fastq-output \\
        ${kreportFlag} \\
        ${args}
    gzip *.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        extract_kraken_reads.py: ${VERSION}
    END_VERSIONS
    """
}
