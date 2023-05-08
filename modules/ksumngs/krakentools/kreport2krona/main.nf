def VERSION = '1.2' // Version information not provided by tool on CLI

process KRAKENTOOLS_KREPORT2KRONA {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::krakentools=1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakentools:1.2--pyh5e36f6f_0':
        'quay.io/biocontainers/krakentools:1.2--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(kreport)

    output:
    tuple val(meta), path("*.krona"), emit: krona
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kreport2krona.py \\
        -r ${kreport} \\
        -o ${prefix}.krona \\
        ${args}

    # Remove ugly 'x__' prefixes for each of the taxonomic levels
    LEVELS=(d k p c o f g s)
    for L in "\${LEVELS[@]}"; do
        sed -i "s/\${L}__//g" ${prefix}.krona
    done

    # Remove underscores that are standing in place of spaces
    sed -i "s/_/ /g" ${prefix}.krona

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kreport2krona.py: ${VERSION}
    END_VERSIONS
    """
}
