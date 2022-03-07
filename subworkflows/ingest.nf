#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CAT_FASTQ } from '../modules/ksumngs/nf-modules/cat/fastq/main.nf'
include { SEQKIT_SPLIT2 } from '../modules/nf-core/modules/seqkit/split2/main.nf'

/// summary: |
///   Take reads from the input folder or a samplesheet and reformat them to be
///   single files with clean names
/// output:
///   - tuple:
///       - type: val(String)
///         description: Sample identifier
///       - type: path
///         description: Reads files
workflow READS_INGEST {
    main:
    versions = Channel.empty()

    // First sanity check: --input must exist
    if (!file(params.input).exists()) {
        log.error "ERROR: file or directory '${params.input}' does not exist!"
        exit 1
    }

    // Sanity check: interleaved reads cannot be single-end
    if (params.interleaved && !params.paired) {
        log.error "ERROR: --interleaved cannot be specified if --paired is false"
        exit 1
    }

    if (file(params.input).isFile()) {
        // --input represents a samplesheet

        // Parse the samplesheet to nf-core sample channel
        Channel
            .of(file("${params.input}"))
            .splitCsv(sep: "\t")
            .filter { !(it[0] ==~ /^#.*/) }
            .map {
                [
                    ['id': it[0], 'single_end': !(params.paired && !params.interleaved), 'strandedness': null],
                    files_transform(it.drop(1))
                ]
            }
            .set { ListedSamples }
    }
    else if (file(params.input).isDirectory()) {
        // --input represents a directory of reads
        if (params.paired && !params.interleaved) {
            // Paired reads can be directly sent to interleaved jail
            Channel
                .fromFilePairs("${params.input}/*{R1,R2,_1,_2}*.{fastq,fq,fastq.gz,fq.gz}")
                .map { [
                        [ 'id': it[0].split('_')[0], 'single_end': false, 'strandedness': null ],
                        it[1]
                        ]
                    }
                .set { ListedSamples }
        }
        else {
            // Reformat single-end/interleaved reads
            Channel
                .fromPath("${params.input}/*.{fastq,fq,fastq.gz,fq.gz}")
                .map{ file -> tuple(file.simpleName, file) }
                .map{ [
                    [ 'id': it[0].split('_')[0], 'single_end': true, 'strandedness': null ],
                    it[1]
                ] }
                .set { ListedSamples }
        }
    }

    // Concatenate the reads together
    CAT_FASTQ(ListedSamples, true)
    CAT_FASTQ.out.reads.set { InterleavedSamples }
    versions = versions.mix(CAT_FASTQ.out.versions)

    if (params.interleaved) {
        // Deinterleave any interleaved reads
        SEQKIT_SPLIT2(InterleavedSamples)
        SEQKIT_SPLIT2.out.reads
            .map { [
                ['id': it[0]['id'], 'single_end': false, 'strandedness': null ],
                it[1]
            ] }
            .set { sample_info }
        versions = versions.mix(SEQKIT_SPLIT2.out.versions)
    }
    else {
        // Transfer reads out of interleaved jail and output them
        InterleavedSamples.set { sample_info }
    }

    emit:
    sample_info
    versions
}

/// summary: |
///   Takes a list of reads files, and returns a new list of file objects
def files_transform(List files) {
    readFiles = [];
    for (int i = 0; i <= files.size(); i+=1) {
        if (files[i]?.trim()) {
            filepath = file(files[i])
            if ( filepath instanceof List ) {
                filepath.each { readFiles.add(it) }
            }
            else {
                readFiles.add(filepath)
            }
        }
    }
    return readFiles
}
