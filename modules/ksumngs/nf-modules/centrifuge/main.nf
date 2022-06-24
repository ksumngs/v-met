process CENTRIFUGE {
    tag "$meta.id"
    label 'process_high'
    label 'process_high_memory'

    conda (params.enable_conda ? 'bioconda::centrifuge=1.0.4_beta' : null)
    
    input:
    tuple val(meta), path(reads)
    path(db)

    output:

    when
    task.ext.when == null || task.ext.when

    script:
    
}