#!/usr/bin/env nextflow

if (params.help) {
    log.info \
    """
NAME
    viral-metagenomics-pipeline - Automated analysis of viral reads in metagenomics samples

SYNOPSIS
    nextflow run millironx/viral-metagenomics-pipeline
        --kraken.db <kraken2 database location>

OPTIONS
    --readsfolder
        The folder containing parired-end Illumina reads in gzipped fastq format. Defaults
        to the current directory

    --threads
        Number of threads to process each sample with. Can't be adjusted on a per-process
        basis. Defaults to 4

    --outfolder
        The place where the final anlysis products will be stored. Defaults to runname_out

    --dev
        Run using fewer inputs and faster process options

    --devinputs
        The number of inputs to take in when using --dev

PROCESS-SPECIFIC OPTIONS
Trimmomatic:
    Please see https://github.com/usadellab/Trimmomatic for full documentation and
    descriptions of the trimming steps.

    --trim-adapters
        Passed to ILLUMINACLIP trimming step. Specifies the path to a fasta file containing
        all the adapters. Valid values are:
            'NexteraPE-PE.fa'
            'TruSeq2-PE.fa'
            'TruSeq3-PE-2.fa'
            'TruSeq3-PE.fa'
        Defaults to 'NexteraPI-PE.fa'. If custom adapters are desired, please see
        Trimmomatic documentation and ensure that the location of the adapter is available
        to the container running Trimmomatic

    --trim-mismatches
        Passed to ILLUMINACLIP trimming step. Specifies the maximum mismatch count which
        will still allow a full adapter match. Defaults to 2

    --trim-pclip
        Passed to ILLUMINACLIP trimming step. Specifies how accurate the match between the
        two reads must be for PE palindrome read alignment. Defaults to 30

    --trim-clip
        Passed to ILLUMINACLIP trimming step. Specifies how accurate the match between any
        adapter sequence must be against a read. Defaults to 10

    --trim-winsize
        Passed to SLIDINGWINDOW trimming step. Specifies the number of bases to average
        across. If used, --trim-winqual must also be specified.

    --trim-winqual
        Passed to the SLIDINGWINDOW trimming step. Specifies the average base quality
        required. If used, --trim-winsize must also be specified.

    --trim-leading
        Passed to the LEADING trimming step. Specifies the minimum quality required to keep
        a base

    --trim-trailing
        Passed to the TRAILING trimming step. Specifies the minimum quality required to keep
        a base

    --trim-crop
        Passed to the CROP trimming step. The number of bases to keep, from the start of
        the read

    --trim-headcrop
        Passed to the HEADCROP trimming step. The number of bases to remove from the start
        of the read

    --trim-minlen
        Passed to the MINLEN trimming step. Specifies the minimum length of reads to
        be kept.
"""
exit 0
}

// Make params persist that need to
RunName = params.runname

// Create an outfolder name if one wasn't provided
if(params.outfolder == "") {
    OutFolder = RunName + "_out"
}
else {
    OutFolder = params.outfolder
}

// Bring in the reads files
Channel
    .fromFilePairs("${params.readsfolder}/*{R1,R2,_1,_2}*.{fastq,fq}.gz")
    .take( params.dev ? params.devinputs : -1 )
    .set{ RawReads }

// First trim, using Trimmomatic
process trim {
    cpus params.threads

    publishDir OutFolder, mode: 'copy'

    input:
    set val(sampleName), file(readsFiles) from RawReads

    output:
    tuple sampleName, file("${sampleName}_trimmomatic_{R1,R2}.fastq.gz") into TrimmedReads
    file("${sampleName}_trimmomatic_{R1,R2}.fastq.gz") into PreKrakenReads

    script:
    // Put together the trimmomatic parameters
    ILLUMINACLIP = "ILLUMINACLIP:/Trimmomatic-0.39/adapters/${params.trimAdapters}:${params.trimMismatches}:${params.trimPclip}:${params.trimClip}"
    SLIDINGWINDOW = ( params.trimWinsize > 0 && params.trimWinqual > 0 ) ? "SLIDINGWINDOW:${params.trimWinsize}:${params.trimWinqual}" : ""
    LEADING = ( params.trimLeading > 0 ) ? "LEADING:${params.trimLeading}" : ""
    TRAILING = ( params.trimTrailing > 0 ) ? "TRAILING:${params.trimTrailing}" : ""
    CROP = ( params.trimCrop > 0 ) ? "CROP:${params.trimCrop}" : ""
    HEADCROP = ( params.trimHeadcrop > 0 ) ? "HEADCROP:${params.trimHeadcrop}" : ""
    MINLEN = ( params.trimMinlen > 0 ) ? "MINLEN:${params.trimMinlen}" : ""
    trimsteps = ILLUMINACLIP + ' ' + SLIDINGWINDOW + ' ' + LEADING + ' ' + TRAILING + ' ' + CROP + ' ' + HEADCROP + ' ' + MINLEN
    """
    trimmomatic PE -threads ${params.threads} \
        ${readsFiles} \
        ${sampleName}_trimmomatic_R1.fastq.gz \
        /dev/null \
        ${sampleName}_trimmomatic_R2.fastq.gz \
        /dev/null \
        ${trimsteps}
    """
}

}
