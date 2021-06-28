#!/usr/bin/env nextflow

/* Parameter Delcarations
    --readsfolder       The folder containing parired-end Illumina reads in
                            gzipped fastq format. Defaults to the current
                            directory
    --threads           Number of threads to process each sample with. Can't be
                            adjusted on a per-process basic. Defaults to 4
    --krakendb          The storage location of the Kraken2 database. Defaults
                            to /kraken2-db
    --blastdb           The storage location of the NCBI BLAST database.
                            Defaults to /blastdb
    --runname           A friendly identifier to describe the samples being
                            analyzed. Defaults to 'viral-metagenomics'
    --outfolder         The place where the final anlysis products will be
                            stored. Defaults to runname_out
*/
params.readsfolder = "."
params.threads = 4
params.outfolder = ""
params.krakendb = "/kraken2-db"
params.blastdb = "/blastdb"
params.runname = "viral-metagenomics"
params.dev = false
params.devinputs = 2

// Make params persist that need to
NumThreads = params.threads
KrakenDb = params.krakendb
RunName = params.runname
BlastDb = params.blastdb
BlastAlgorithms = ['blastn', 'blastx']

// Create an outfolder name if one wasn't provided
if(params.outfolder == "") {
    OutFolder = RunName + "_out"
}
else {
    OutFolder = params.outfolder
}

// Bring in the reads files
Channel
    .fromFilePairs("${params.readsfolder}/*{R1,R2,_1,_2}*.{fastq,fq}.{gz,bz,bz2}")
    .take( params.dev ? params.devinputs : -1 )
    .set{ RawReads }

// First trim, using Trimmomatic
process trimmomatic {
    input:
    set val(sampleName), file(readsFiles) from RawReads

    output:
    tuple sampleName, file("${sampleName}_trimmomatic_{R1,R2}.fastq.gz") into OnceTrimmedReads

    script:
    """
    # Very generic trimmomatic settings, except that unpaired reads are discarded
    trimmomatic PE -threads ${NumThreads} \
        ${readsFiles} \
        ${sampleName}_trimmomatic_R1.fastq.gz \
        /dev/null \
        ${sampleName}_trimmomatic_R2.fastq.gz \
        /dev/null \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Secord trim, using SeqPurge
process seqpurge {
    input:
    set val(sampleName), file(readsFiles) from OnceTrimmedReads

    // These trimmed reads start a branching path:
    // First branch: get classified by Kraken and produce visuals
    // Secord branch: based on classification, go to de novo assembly
    output:
    tuple sampleName, file("${sampleName}_seqpurge_{R1,R2}.fastq.gz") into TwiceTrimmedReads
    file("${sampleName}_seqpurge_{R1,R2}.fastq.gz") into KrakenInputReads

    script:
    """
    SeqPurge -threads ${NumThreads} \
        -in1  ${readsFiles[0]} \
        -in2  ${readsFiles[1]} \
        -out1 ${sampleName}_seqpurge_R1.fastq.gz \
        -out2 ${sampleName}_seqpurge_R2.fastq.gz
    """
}

// Classify reads using Kraken
process kraken {
    input:
    set val(sampleName), file(readsFiles) from TwiceTrimmedReads

    output:
    tuple sampleName, file("${sampleName}.kraken"), file("${sampleName}.krpt") into KrakenFile
    tuple sampleName, file("${sampleName}.krpt") into KrakenVisuals

    script:
    quickflag = params.dev ? '--quick' : ''
    """
    kraken2 --db ${KrakenDb} --threads ${NumThreads} --paired ${quickflag} \
        --report "${sampleName}.krpt" \
        --output "${sampleName}.kraken" \
        ${readsFiles}
    """
}

// Pull the viral reads and any unclassified reads from the original reads
// files for futher downstream processing using KrakenTools
process filterreads {
    input:
    file(readsFiles) from KrakenInputReads
    set val(sampleName), file(krakenFile), file(krakenReport) from KrakenFile

    output:
    tuple sampleName, file("${sampleName}_filtered_{R1,R2}.fastq.gz") into FilteredReads
    tuple sampleName, file("${sampleName}_filtered_{R1,R2}.fastq.gz") into FilteredReads2
    tuple sampleName, file("${sampleName}_filtered_{R1,R2}.fastq.gz") into FilteredReads3

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

// Convert the kraken reports to krona input files using KrakenTools then
// prettify them using unix tools
process kraken2krona {
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
    // This is a final step: publish it
    publishDir OutFolder, mode: 'move'

    input:
    file '*' from KronaText.collect()

    output:
    file("${RunName}.html")

    // Name the file using the run name, and name the top level taxa as 'root'
    // consistent with kraken
    script:
    """
    ktImportText * -o ${RunName}.html -n root
    """
}

// Assemble into contigs using Ray
process ray {
    input:
    set val(sampleName), file(readsFiles) from FilteredReads

    output:
    tuple val(sampleName), val(assembler), 'RayOutput/Contigs.fasta' into RayContigs

    script:
    assembler = 'ray'
    """
    mpiexec -n ${NumThreads} Ray -p ${readsFiles}
    """

}

// Assemble into contigs using iva
process iva {
    input:
    set val(sampleName), file(readsFiles) from FilteredReads2

    output:
    tuple val(sampleName), val(assembler), 'contigs.fasta' into IVAContigs

    script:
    assembler = 'iva'
    """
    iva -t ${NumThreads} -f ${readsFiles[0]} -r ${readsFiles[1]} out
    mv out/contigs.fasta .
    """
}

// Assemble into contigs using a5
process a5 {
    input:
    set val(sampleName), file(readsFiles) from FilteredReads3

    output:
    tuple val(sampleName), val(assembler), "${sampleName}.contigs.fasta" into A5Contigs

    script:
    assembler = 'a5'
    """
    a5_pipeline.pl --threads ${NumThreads} ${readsFiles} ${sampleName}
    """
}

// Blast contigs
process blast {
    // This is a final step: publish it
    publishDir OutFolder, mode: 'copy'

    // Blast needs to happen on all contigs from all assemblers, and both
    // blastn and blastx needs to be applied to all contigs
    input:
    set val(sampleName), val(assembler), file(readsFiles) from RayContigs.concat(IVAContigs, A5Contigs)
    each program from BlastAlgorithms

    output:
    file("${RunName}_${sampleName}_${assembler}.${program}.tsv")

    script:
    // Pick the faster algorithm if this is a development cycle, otherwise
    // the titular program is also the name of the algorithm
    algorithm = program
    if ( params.dev ) {
        algorithm = ( program == 'blastn' ) ? 'megablast' : 'blastx-fast'
    }

    // Switch which database to read from
    dbExtension = (program == 'blastn') ? 'nt' : 'nr'

    // Take only first hits if this is a development cycle
    maxhsps = params.dev ? 1 : 10
    nalignments = params.dev ? 1 : 5

    // Squash the filename into a single variable
    outFile = "${RunName}_${sampleName}_${assembler}.${program}.tsv"
    """
    echo "Sequence ID\tDescription\tGI\tTaxonomy ID\tScientific Name\tCommon Name\tRaw score\tBit score\tQuery Coverage\tE value\tPercent identical\tSubject length\tAlignment length\tAccession\tMismatches\tGap openings\tStart of alignment in query\tEnd of alignment in query\tStart of alignment in subject\tEnd of alignment in subject" > ${outFile}
    ${program} -query ${readsFiles} \
        -db ${BlastDb}/${dbExtension} \
        -max_hsps ${maxhsps} \
        -num_alignments ${nalignments} \
        -outfmt "6 qseqid stitle sgi staxid ssciname scomname score bitscore qcovs evalue pident length slen saccver mismatch gapopen qstart qend sstart send" \
        -evalue 1e-5 \
        -num_threads ${NumThreads} \
        -task ${algorithm} >> ${outFile}
    """
}
