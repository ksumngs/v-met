// Adapted from nf-core's Kraken2/Kraken2 module by Jose Espinosa-Carrasco and
// Harshil Patel
// (https://github.com/nf-core/modules/blob/master/modules/kraken2/kraken2)
// Used under MIT license
// This module is modified from the original to include a high-memory label, and to
// emit the kraken output as a gzipped text file
process KRAKEN2 {
    tag "$meta.id"
    label 'process_high'
    label 'process_high_memory'

    conda (params.enable_conda ? 'bioconda::kraken2=2.1.2 conda-forge::pigz=2.6' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' :
        'quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' }"

    input:
    tuple val(meta), path(reads)
    path(db)

    output:
    tuple val(meta), path("*classified*")  , emit: classified
    tuple val(meta), path("*unclassified*"), emit: unclassified
    tuple val(meta), path("*.kraken.gz")   , emit: kraken
    tuple val(meta), path("*.kreport")     , emit: kreport
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pairedFlag = meta.single_end ? '' : '--paired'
    def classifiedFlag = meta.single_end ? "${prefix}_classified.fastq" : "${prefix}_classified#.fastq"
    def unclassifiedFlag = meta.single_end ? "${prefix}_unclassified.fastq" : "${prefix}_unclassified#.fastq"
    """
    kraken2 \\
            --db ${db} \\
            --threads ${task.cpus} \\
            --classified-out ${classifiedFlag} \\
            --unclassified-out ${unclassifiedFlag} \\
            --report ${prefix}.kreport \\
            ${pairedFlag} \\
            ${args} \\
            ${reads} \\
        | gzip \\
        > ${prefix}.kraken.gz

    pigz -p${task.cpus} *.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
