#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTQC } from '../modules/nf-core/modules/fastqc/main.nf'
include { NANOSTAT } from '../modules/ksumngs/nf-modules/nanostat/main.nf'
include { SEQTK_MERGEPE } from '../modules/nf-core/modules/seqtk/mergepe/main.nf'

/// summary: |
///   Perform context-sensitive QC on fastq reads
workflow QC {
    take:
    reads

    main:
    versions = Channel.empty()

    if (params.platform == 'illumina') {
        SEQTK_MERGEPE(reads)
        FASTQC(SEQTK_MERGEPE.out.reads)
        FASTQC.out.zip.set{ report }

        versions = versions.mix(SEQTK_MERGEPE.out.versions)
        versions = versions.mix(FASTQC.out.versions)
    }
    else if (params.platform == 'nanopore') {
        NANOSTAT(reads)
        NANOSTAT.out.log.set{ report }

        versions = versions.mix(NANOSTAT.out.versions)
    }

    emit:
    report
    versions
}
