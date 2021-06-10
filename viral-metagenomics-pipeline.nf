#!/usr/bin/env nextflow

params.readsfolder
params.threads

NumThreads = params.threads

Channel
    .fromFilePairs("${params.readsfolder}/*{R1,R2,_1,_2}*.{fastq,fq,fasta,fa}.{gz,bz,bz2}")
    .set{ RawReads }

process trimmomatic {
    input:
        set val(sampleName), file(readsFiles) from RawReads

    output:
        tuple sampleName, file("${sampleName}_trimmomatic_{R1,R2}.fastq.gz") into TrimmomaticReads

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
        set val(sampleName), file(readsFiles) from TrimmomaticReads

    output:
        tuple sampleName, file("${sampleName}_seqpurge_{R1,R2}.fastq.gz") into SeqPurgeReads

    script:
    """
        SeqPurge -threads ${NumThreads} \
            -in1  ${readsFiles[0]} \
            -in2  ${readsFiles[1]} \
            -out1 ${sampleName}_seqpurge_R1.fastq.gz \
            -out2 ${sampleName}_seqpurge_R2.fastq.gz
    """
}
