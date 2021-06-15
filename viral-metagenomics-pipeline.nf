#!/usr/bin/env nextflow

params.readsfolder = "."
params.threads = 4
params.krakendb = "/kraken2-db"

NumThreads = params.threads
KrakenDb = params.krakendb

Channel
    .fromFilePairs("${params.readsfolder}/*{R1,R2,_1,_2}*.{fastq,fq,fasta,fa}.{gz,bz,bz2}")
    .set{ RawReads }

process trimmomatic {
    input:
        set val(sampleName), file(readsFiles) from RawReads

    output:
        tuple sampleName, file("${sampleName}_trimmomatic_{R1,R2}.fastq.gz") into OnceTrimmedReads

        script:
        """
            trimmomatic PE -threads ${NumThreads} \
                ${readsFiles} \
                ${sampleName}_trimmomatic_R1.fastq.gz \
                /dev/null \
                ${sampleName}_trimmomatic_R2.fastq.gz \
                /dev/null \
                ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """
}

process seqpurge {
    conda 'bioconda::ngs-bits'

    input:
        set val(sampleName), file(readsFiles) from OnceTrimmedReads

    output:
        tuple sampleName, file("${sampleName}_seqpurge_{R1,R2}.fastq.gz") into TwiceTrimmedReads
        tuple sampleName, file("${sampleName}_seqpurge_{R1,R2}.fastq.gz") into KrakenInputReads

    script:
    """
        SeqPurge -threads ${NumThreads} \
            -in1  ${readsFiles[0]} \
            -in2  ${readsFiles[1]} \
            -out1 ${sampleName}_seqpurge_R1.fastq.gz \
            -out2 ${sampleName}_seqpurge_R2.fastq.gz
    """
}

process kraken {
    input:
        set val(sampleName), file(readsFiles) from TwiceTrimmedReads

    output:
        tuple sampleName, file("${sampleName}.kraken"), file("${sampleName}.krpt") into KrakenFile
        tuple sampleName, file("${sampleName}.krpt") into KrakenVisuals

    script:
    """
        kraken2 --db ${KrakenDb} --threads ${NumTehreads} --paired \
            --use-names \
            --report "${sampleName}.krpt" \
            --output "${sampleName}.kraken" \
            ${readsFiles}
    """

}

process filterreads {
    conda 'bioconda::krakentools'

    input:
        set val(devnull), file(readsFiles) from KrakenInputReads
        set val(sampleName), file(krakenFile), file(krakenReport) from KrakenFile

    output:
        tuple sampleName, file("${sampleName}_filtered_{R1,R2}.fastq.gz") into FilteredReads

    // Although I haven't seen it documented anywhere, 0 is unclassified reads
    // and 10239 is viral reads
    script:
    """
        extract_kraken_reads.py -k ${krakenFile} \
        -s1 ${readsFiles[0]} -s2 ${readsFiles[1]} \
        -r ${krakenReport} \
        -t 0 10239 --include-children \
        --fastq-output \
        -o ${sampleName}_filtered_R1.fastq -o2 ${sampleName}_filtered_R2.fastq
        gzip ${sampleName}_filtered_{R1,R2}.fastq
    """

}

process kraken2krona {
    conda 'bioconda::krakentools'

    input:
        set val(sampleName), file(krakenReport) from KrakenVisuals

    output:
        tuple val(sampleName), file("${sampleName}.krona") into KronaText

    shell:
    '''
        kreport2krona.py -r !{krakenReport} -o !{sampleName}.krona
        LEVELS=(d k p c o f g s)
        for L in $LEVELS; do
            sed -i "s/${L}__//" !{sampleName}.krona
        done
    '''
}

process krona {
    conda 'bioconda::krona'

    input:
        set val(sampleName), file(kronaText) from KronaText

    output:
        tuple val(sampleName), file("${sampleName}.html") into KronaWebPage

    script:
    """
        ktImportText ${kronaText} -o ${sampleName}.html
    """
}

/*
process ray {
    input:
        set val(sampleName), file(readsFiles) from FilteredReads

    output:
        tuple val(sampleName), file("RayOutput/Contigs.fasta") into RayContigs
        tuple val(sampleName), file("RayOutput/Scaffolds.fasta") into RayScaffolds

    script:
    """
        mpiexec -n ${NumThreads} Ray -p ${readsFiles}
    """

}
*/
