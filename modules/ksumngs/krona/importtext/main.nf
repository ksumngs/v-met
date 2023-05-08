def VERSION = '2.8.1' // Version information not provided by tool on CLI

process KRONA_IMPORTTEXT {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::krona=2.8.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krona:2.8.1--pl5321hdfd78af_1':
        'quay.io/biocontainers/krona:2.8.1--pl5321hdfd78af_1' }"

    input:
    path '*'

    output:
    path "krona.html", emit: html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ktImportText * \\
        -o krona.html \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        KronaTools: ${VERSION}
    END_VERSIONS
    """
}
