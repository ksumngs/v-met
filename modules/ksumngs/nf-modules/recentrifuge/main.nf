process RECENTRIFUGE {
    tag "$meta.id"
    label 'process_high'
    label 'process_high_memory'

    conda (params.enable_conda ? 'bioconda::recentrifuge=1.9.0' : null)
    
    input:

    output:

    when
    task.ext.when == null || task.ext.when

    script:
    
}