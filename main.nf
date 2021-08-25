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

    --runname
        A friendly identifier to describe the samples being analyzed. Defaults to
        'viral-metagenomics'

    --outfolder
        The place where the final anlysis products will be stored. Defaults to runname_out

    --dev
        Run using fewer inputs and faster process options

    --devinputs
        The number of inputs to take in when using --dev

    --kmer-length
        The length of kmers when assembling. Defaults to 35

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

Kraken:
    See https://github.com/DerrickWood/kraken2/wiki/Manual for full documentation of
    Kraken 2's available options
    --kraken-db
        Path to Kraken 2 database. REQUIRED

BLAST:
    --blast-db
        Path to blast databases. REQUIRED. It is also recommended to place this value in the
        BLASTDB environment variable.
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

// Classify reads using Kraken
process kraken {
    cpus params.threads

    input:
    set val(sampleName), file(readsFiles) from TrimmedReads

    output:
    tuple sampleName, file("${sampleName}.kraken"), file("${sampleName}.krpt") into KrakenFile
    tuple sampleName, file("${sampleName}.krpt") into KrakenVisuals

    script:
    quickflag = params.dev ? '--quick' : ''
    """
    kraken2 --db ${params.krakenDb} --threads ${params.threads} --paired ${quickflag} \
        --report "${sampleName}.krpt" \
        --output "${sampleName}.kraken" \
        ${readsFiles}
    """
}

// Convert the kraken reports to krona input files using KrakenTools then
// prettify them using unix tools
process kraken2krona {
    cpus 1

    input:
    set val(sampleName), file(krakenReport) from KrakenVisuals

    output:
    file("${sampleName}.krona") into KronaText

    shell:
    '''
    #!/bin/bash
    # Using bash-specific loop syntax here, so shebang is required

    # Convert the report using KrakenTools
    kreport2krona.py -r !{krakenReport} -o !{sampleName}.krona

    # KrakenTools creates ugly x__ prefixes for each of the taxonomic levels:
    # let's remove each of those
    LEVELS=(d k p c o f g s)
    for L in "${LEVELS[@]}"; do
        sed -i "s/${L}__//g" !{sampleName}.krona
    done

    # Also remove underscores that are standing in place of spaces
    sed -i "s/_/ /g" !{sampleName}.krona
    '''
}

// Collect all of the krona input files and convert them to a single graphical
// webpage to ship to the user
process krona {
    cpus 1

    // This is a final step: publish it
    publishDir OutFolder, mode: 'copy'

    input:
    file '*' from KronaText.collect()

    output:
    file("krona.html")

    // Name the file using the run name, and name the top level taxa as 'root'
    // consistent with kraken
    script:
    """
    ktImportText * -o krona.html -n root
    """
}

// Pull the viral reads and any unclassified reads from the original reads
// files for futher downstream processing using KrakenTools
process filterreads {
    cpus 1

    input:
    file(readsFiles) from PreKrakenReads
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

// Convert the viral and unclassified reads to fasta for blast
process convert2fasta {
    cpus 1

    input:
    set val(sampleName), file(readsFiles) from FilteredReads

    output:
    tuple val(sampleName), file("${sampleName}.fasta") into FilteredFastas

    script:
    """
    seqtk mergepe ${readsFiles} | seqtk seq -a - > ${sampleName}.fasta
    """
}

// Blast contigs
process blast {
    cpus params.threads

    // This is a final step: publish it
    publishDir OutFolder, mode: 'copy'

    // Blast needs to happen on all contigs from all assemblers, and both
    // blastn and blastx needs to be applied to all contigs
    input:
    set val(sampleName), file(readsFiles) from FilteredFastas

    output:
    file("${sampleName}.blast.tsv")

    script:
    // Separate parameters from script
    max_hsps       = 10
    num_alignments = 5
    outfmt         = '"6 qseqid stitle sgi staxid ssciname scomname score bitscore qcovs evalue pident length slen saccver mismatch gapopen qstart qend sstart send"'
    evalue         = 1e-5

    // Pick the faster algorithm if this is a development cycle, otherwise
    // the titular program is also the name of the algorithm
    if ( params.dev ) {
        max_hsps = 1
        num_alignments = 1
        evalue = 1e-50
    }

    // Squash the filename into a single variable
    outFile = "${sampleName}.blast.tsv"
    """
    echo "Sequence ID\tDescription\tGI\tTaxonomy ID\tScientific Name\tCommon Name\tRaw score\tBit score\tQuery Coverage\tE value\tPercent identical\tSubject length\tAlignment length\tAccession\tMismatches\tGap openings\tStart of alignment in query\tEnd of alignment in query\tStart of alignment in subject\tEnd of alignment in subject" > ${outFile}
    blastn -query ${readsFiles} \
        -db ${params.blastDb}/nt \
        -max_hsps ${max_hsps} \
        -num_alignments ${num_alignments} \
        -outfmt ${outfmt} \
        -evalue ${evalue} \
        -num_threads ${params.threads} >> ${outFile}
    """
}
