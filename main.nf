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
    conda 'bioconda::ngs-bits'

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
    if ( params.dev )
    """
    kraken2 --db ${KrakenDb} --threads ${NumThreads} --paired --quick \
        --report "${sampleName}.krpt" \
        --output "${sampleName}.kraken" \
        ${readsFiles}
    """
    else
    """
    kraken2 --db ${KrakenDb} --threads ${NumThreads} --paired \
        --report "${sampleName}.krpt" \
        --output "${sampleName}.kraken" \
        ${readsFiles}
    """

}

// Pull the viral reads and any unclassified reads from the original reads
// files for futher downstream processing using KrakenTools
process filterreads {
    conda 'bioconda::krakentools'

    input:
        file(readsFiles) from KrakenInputReads
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

// Convert the kraken reports to krona input files using KrakenTools then
// prettify them using unix tools
process kraken2krona {
    conda 'bioconda::krakentools'

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
    conda 'bioconda::krona'

    // This is a final step: publish it
    publishDir OutFolder, mode: 'move'

    input:
        file '*' from KronaText.collect()

    output:
        file("${RunName}.html") into KronaWebPage

    // Name the file using the run name, and name the top level taxa as 'root'
    // consistent with kraken
    script:
    """
    ktImportText * -o ${RunName}.html -n root
    """
}

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

process blast {
    // This is a final step: publish it
    publishDir OutFolder, mode: 'copy'

    input:
        set val(sampleName), file(readsFiles) from RayContigs

    output:
        file("${RunName}_${sampleName}.blast.tsv")

    script:
    if ( params.dev )
    """
    echo "Sequence ID\tDescription\tGI\tTaxonomy ID\tScientific Name\tCommon Name\tRaw score\tBit score\tQuery Coverage\tE value\tPercent identical\tSubject length\tAlignment length\tAccession\tMismatches\tGap openings\tStart of alignment in query\tEnd of alignment in query\tStart of alignment in subject\tEnd of alignment in subject" > ${RunName}_${sampleName}.blast.tsv
    blastn -query ${readsFiles} \
            -db ${BlastDb}/nt \
            -max_hsps 1 \
            -num_alignments 1 \
            -outfmt "6 qseqid stitle sgi staxid ssciname scomname score bitscore qcovs evalue pident length slen saccver mismatch gapopen qstart qend sstart send" \
            -evalue 1e-5 \
            -num_threads ${NumThreads} \
            -task megablast >> ${RunName}_${sampleName}.blast.tsv
    """
    else
    """
    echo "Sequence ID\tDescription\tGI\tTaxonomy ID\tScientific Name\tCommon Name\tRaw score\tBit score\tQuery Coverage\tE value\tPercent identical\tSubject length\tAlignment length\tAccession\tMismatches\tGap openings\tStart of alignment in query\tEnd of alignment in query\tStart of alignment in subject\tEnd of alignment in subject" > ${RunName}_${sampleName}.blast.tsv
    blastn -query ${readsFiles} \
            -db ${BlastDb}/nt \
            -max_hsps 10 \
            -num_alignments 5 \
            -outfmt "6 qseqid stitle sgi staxid ssciname scomname score bitscore qcovs evalue pident length slen saccver mismatch gapopen qstart qend sstart send" \
            -evalue 1e-5 \
            -num_threads ${NumThreads} \
            -task blastn >> ${RunName}_${sampleName}.blast.tsv
    """

}
