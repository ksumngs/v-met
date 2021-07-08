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
    .fromPath("${params.readsfolder}/*.{fastq,fq}.gz")
    .take( params.dev ? params.devinputs : -1 )
    .into{ RawReads; UnclassifiedReads }

// Classify reads using Kraken
process kraken {
    cpus params.threads

    input:
    set file(readsFile) from RawReads

    output:
    tuple sampleName, file("${sampleName}.kraken"), file("${sampleName}.krpt") into KrakenFile
    tuple sampleName, file("${sampleName}.krpt") into KrakenVisuals

    script:
    quickflag = params.dev ? '--quick' : ''
    fileName = readsFile.tokenize('.')
    sampleName = fileName[0]
    """
    kraken2 --db ${params.krakenDb} --threads ${params.threads} ${quickflag} \
        --report "${sampleName}.krpt" \
        --output "${sampleName}.kraken" \
        ${readsFile}
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
    file(readsFile) from UnclassifiedReads
    set val(sampleName), file(krakenFile), file(krakenReport) from KrakenFile

    output:
    tuple sampleName, file("${sampleName}_filtered.fastq.gz") into ReadsForRay
    tuple sampleName, file("${sampleName}_filtered.fastq.gz") into ReadsForMetaVelvet
    tuple sampleName, file("${sampleName}_filtered.fastq.gz") into ReadsForAbyss
    file("${sampleName}_filtered_.fastq.gz") into CompressedReadsForRemapping

    // Although I haven't seen it documented anywhere, 0 is unclassified reads
    // and 10239 is viral reads
    script:
    """
    extract_kraken_reads.py -k ${krakenFile} \
        -s ${readsFile} \
        -r ${krakenReport} \
        -t 0 10239 --include-children \
        --fastq-output \
        -o ${sampleName}_filtered.fastq
    gzip ${sampleName}_filtered.fastq
    """

}

// Assemble into contigs using Ray
process ray {
    cpus params.threads

    input:
    set val(sampleName), file(readsFile) from ReadsForRay

    output:
    tuple val(sampleName), val(assembler), 'Contigs.fasta' into RayContigsForBlast
    tuple val(sampleName), val(assembler), file('Contigs.fasta'), file(readsFile) into RayContigsForRemapping

    script:
    // Export the assembler for future combined steps
    assembler = 'ray'
    """
    mpiexec -n ${params.threads} Ray -k ${params.kmerLength} -s ${readsFile}
    mv RayOutput/Contigs.fasta .
    """
}

// Assemble using MetaVelvet
process metavelvet {
    cpus params.threads

    input:
    set val(sampleName), file(readsFile) from ReadsForMetaVelvet

    output:
    tuple val(sampleName), val(assembler), 'meta-velvetg.contigs.fa' into MetaVelvetContigsForBlast
    tuple val(sampleName), val(assembler), file('meta-velvetg.contigs.fa'), file(readsFile) into MetaVelvetContigsForRemapping

    script:
    assembler = 'metavelvet'
    """
    export OMP_NUM_THREADS=${params.threads}
    export OMP_THREAD_LIMIT=${params.threads}
    velveth out ${params.kmerLength} -fastq.gz -long ${readsFile}
    velvetg out -exp_cov auto
    meta-velvetg out
    mv out/meta-velvetg.contigs.fa .
    """
}

// Assemble using Abyss
process abyss {
    cpus params.threads

    input:
    set val(sampleName), file(readsFile) from ReadsForAbyss

    output:
    tuple val(sampleName), val(assembler), 'contigs.fa' into AbyssContigsForBlast
    tuple val(sampleName), val(assembler), file('contigs.fa'), file(readsFile) into AbyssContigsForRemapping

    script:
    assembler = 'abyss'
    """
    abyss-pe np=${params.threads} name=${sampleName} k=21 se="${readsFile}"
    cp ${sampleName}-contigs.fa contigs.fa
    """
}

// Remap contigs using BWA
process bwa {
    cpus params.threads

    input:
    set val(sampleName), val(assembler), file(contigs), file(readsFile) from RayContigsForRemapping.concat(MetaVelvetContigsForRemapping, AbyssContigsForRemapping)

    output:
    tuple val(sampleName), val(assembler), file(contigs), file("${sampleName}_${assembler}.sam") into RemappedReads

    script:
    """
    cp ${readsFile} read.fastq.gz
    gunzip read.fastq.gz
    bwa index ${contigs}
    bwa aln -t ${params.threads} ${contigs} read.fastq > ${sampleName}.sai
    bwa samse ${contigs} \
        ${sampleName}.sai ${readsFile} > ${sampleName}_${assembler}.sam
    """
}

// Sort and compress the sam files for visualization
process sortsam {
    cpus 1

    input:
    set val(sampleName), val(assembler), file(contigs), file(samfile) from RemappedReads

    output:
    tuple file("${sampleName}_${assembler}.contigs.fasta"), file("${sampleName}_${assembler}.contigs.fasta.fai"), file("${sampleName}_${assembler}.bam"), file("${sampleName}_${assembler}.bam.bai") into Assemblies

    script:
    """
    # Rename the contigs file to a consistent format
    mv ${contigs} ${sampleName}_${assembler}.contigs.fasta

    # Create a contigs indes
    samtools faidx ${sampleName}_${assembler}.contigs.fasta

    # Convert and sort the sam file
    samtools view -S -b ${samfile} > sample.bam
    samtools sort sample.bam -o ${sampleName}_${assembler}.bam

    # Index the sorted bam file
    samtools index ${sampleName}_${assembler}.bam
    """
}

// Create a viewer of all the assembly files
process assemblyview {
    cpus 1

    publishDir OutFolder, mode: 'copy'

    input:
    file '*' from Assemblies.collect()

    output:
    file 'index.html'
    file 'index.js'
    file 'package.json'
    file 'data/*'

    script:
    """
    mkdir data
    mv *.contigs.fasta *.contigs.fasta.fai *.bam *.bam.bai data
    git clone https://github.com/MillironX/igv-bundler.git igv-bundler
    mv igv-bundler/{index.html,index.js,package.json} .
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
    set val(sampleName), val(assembler), file(readsFiles) from RayContigsForBlast.concat(MetaVelvetContigsForBlast, AbyssContigsForBlast)

    output:
    file("${RunName}_${sampleName}_${assembler}.blast.tsv")

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
    outFile = "${RunName}_${sampleName}_${assembler}.blast.tsv"
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
