#!/usr/bin/env nextflow

/*
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

    --blastdb
        The storage location of the NCBI BLAST database. Defaults to /blastdb

    --runname
        A friendly identifier to describe the samples being analyzed. Defaults to
        'viral-metagenomics'

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

    --trimmomatic.fastaWithAdapters
        Passed to ILLUMINACLIP trimming step. Specifies the path to a fasta file containing
        all the adapters. Valid values are:
            'NexteraPE-PE.fa'
            'TruSeq2-PE.fa'
            'TruSeq3-PE-2.fa'
            'TruSeq3-PE.fa'
        Defaults to 'NexteraPI-PE.fa'. If custom adapters are desired, please see
        Trimmomatic documentation and ensure that the location of the adapter is available
        to the container running Trimmomatic

    --trimmomatic.seedMismatches
        Passed to ILLUMINACLIP trimming step. Specifies the maximum mismatch count which
        will still allow a full adapter match. Defaults to 2

    --trimmomatic.palindromeClipThreshold
        Passed to ILLUMINACLIP trimming step. Specifies how accurate the match between the
        two reads must be for PE palindrome read alignment. Defaults to 30

    --trimmomatic.simpleClipThreshold
        Passed to ILLUMINACLIP trimming step. Specifies how accurate the match between any
        adapter sequence must be against a read. Defaults to 10

    --trimmomatic.windowSize
        Passed to SLIDINGWINDOW trimming step. Specifies the number of bases to average
        across. If used, --trimmomatic.requiredQuality must also be specified.

    --trimmomatic.requiredQuality
        Passed to the SLIDINGWINDOW trimming step. Specifies the average base quality
        required. If used, --trimmomatic.windowSize must also be specified.

    --trimmomatic.leading
        Passed to the LEADING trimming step. Specifies the minimum quality required to keep
        a base

    --trimmomatic.trailing
        Passed to the TRAILING trimming step. Specifies the minimum quality required to keep
        a base

    --trimmomatic.crop
        Passed to the CROP trimming step. The number of bases to keep, from the start of
        the read

    --trimmomatic.headcrop
        Passed to the HEADCROP trimming step. The number of bases to remove from the start
        of the read

    --trimmomatic.minlen
        Passed to the MINLEN trimming step. Specifies the minimum length of reads to
        be kept.

SeqPurge:
    See https://github.com/imgag/ngs-bits/blob/master/doc/tools/SeqPurge.md for full
    parameter descriptions.
    --seqpurge.a1
        Forward adapter sequence (at least 15 bases). Defaults
        to 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'

    --seqpurge.a2
        Reverse adapter sequence (at least 15 bases). Defaults
        to 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'

    --seqpurge.match_perc
        Minimum percentage of matching bases for sequence/adapter matches. Defaults to 80

    --seqpurge.mep
        Maximum error probability of insert and adapter matches. Defaults to 1e-6

    --seqpurge.qcut
        Quality trimming cutoff for trimming from the end of reads using a sliding window
        approach. Set to 0 to disable. Defaults to 15

    --seqpurge.qwin
        Quality trimming window size. Defaults to 5

    --seqpurge.qoff
        Quality trimming FASTQ score offset. Defaults to 33

    --seqpurge.ncut
        Number of subsequent Ns to trimmed using a sliding window approach from the front of
        reads. Set to 0 to disable. Defaults to 7

    --seqpurge.min_len
        Minimum read length after adapter trimming. Shorter reads are discarded. Defaults
        to 30

Kraken:
    See https://github.com/DerrickWood/kraken2/wiki/Manual for full documentation of
    Kraken 2's available options
    --kraken.db
        Path to Kraken 2 database. REQUIRED

Ray:
    See https://github.com/sebhtml/ray/blob/master/MANUAL_PAGE.txt for full documentation of
    Ray's available options
    --ray.kmerlength
        Selects the length of k-mers. It must be odd because reverse-complement vertices are
        stored together. Larger k-mers utilise more memory. Defaults to 21

    --ray.minimumSeedLength
        Changes the minimum seed length. Defaults to 100

    --ray.minimumContigLength
        Changes the minimum contig length. Defaults to 100

    --ray.maximumSeedCoverageDepth
        Ignores any seed with a coverage depth above this threshold. Defaults to 4294967295

    --ray.minimumSeedCoverageDepth
        Sets the minimum seed coverage depth. Any path with a coverage depth lower than this
        will be discarded. Defaults to 0

BWA:
    See http://bio-bwa.sourceforge.net/bwa.shtml for full documentation of BWA's available
    options. Note that these options only apply to the alignment step.
    --bwa.n
        Maximum edit distance if the value is INT, or the fraction of missing alignments
        given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit
        distance is automatically chosen for different read lengths. Defaults to 0.04

    --bwa.o
        Maximum number of gap opens. Defaults to 1

    --bwa.e
        Maximum number of gap extensions, -1 for k-difference mode (disallowing long gaps).
        Defaults to -1

    --bwa.d
        Disallow a long deletion within so many bp towards the 3'-end. Defaults to 16

    --bwa.i
        Dissalow an indel within so many bp towards the ends. Defaults to 5

    --bwa.k
        Maximum edit distance in the seed

    --bwa.M
        Mismatch penalty. Defaults to 3

    --bwa.O
        Gap open penalty. Defaults to 4

    --bwa.E
        Gap extension penalty. Defaults to 4

BLAST:
    --blast.db
        Path to blast databases. REQUIRED. It is also recommended to place this value in the
        BLASTDB environment variable.

    --blast.max_hsps
        Max hits from blast. Defaults to 10

    --blast.num_alignments
        Max alignments from blast. Defaults to 5

    --blast.outfmt
        Output format of blast. Defaults to
        '"6 qseqid stitle sgi staxid ssciname scomname score bitscore qcovs evalue pident length slen saccver mismatch gapopen qstart qend sstart send"'
        It is highly recommended to not alter this, as a heading will be generated based on
        the default value

    --blast.evalue
        The maximum e-value allowed. Defaults to 1e-5
*/

params.readsfolder = "."
params.threads = 4
params.outfolder = ""
params.runname = "viral-metagenomics"
params.dev = false
params.devinputs = 2

// Make params persist that need to
NumThreads = params.threads
RunName = params.runname
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
    .fromFilePairs("${params.readsfolder}/*{R1,R2,_1,_2}*.{fastq,fq}.gz")
    .take( params.dev ? params.devinputs : -1 )
    .set{ RawReads }

// First trim, using Trimmomatic
process trim {
    input:
    set val(sampleName), file(readsFiles) from RawReads

    output:
    tuple sampleName, file("${sampleName}_trimmomatic_{R1,R2}.fastq.gz") into TrimmedReads

    script:
    // Specifying Trimmomatic paramaters on the command line can sometimes wipe
    // out the IlluminaClip settings, so restore them if absent
    if ( params.trimmomatic.fastaWithAdapters == null ) {
        params.trimmomatic.fastaWithAdapters = 'NexteraPE-PE.fa'
        params.trimmomatic.seedMismatches = 2
        params.trimmomatic.palindromeClipThreshold = 30
        params.trimmomatic.simpleClipThreshold = 10
    }

    // Put together the trimmomatic parameters
    ILLUMINACLIP = "ILLUMINACLIP:${params.trimmomatic.fastaWithAdapters}:${params.trimmomatic.seedMismatches}:${params.trimmomatic.palindromeClipThreshold}:${params.trimmomatic.simpleClipThreshold}"
    SLIDINGWINDOW = ( params.trimmomatic.windowSize > 0 && params.trimmomatic.requiredQuality > 0 ) ? "SLIDINGWINDOW:${params.trimmomatic.windowSize}:${params.trimmomatic.requiredQuality}" : ""
    LEADING = ( params.trimmomatic.leading > 0 ) ? "LEADING:${params.trimmomatic.leading}" : ""
    TRAILING = ( params.trimmomatic.trailing > 0 ) ? "TRAILING:${params.trimmomatic.trailing}" : ""
    CROP = ( params.trimmomatic.crop > 0 ) ? "CROP:${params.trimmomatic.crop}" : ""
    HEADCROP = ( params.trimmomatic.headcrop > 0 ) ? "HEADCROP:${params.trimmomatic.headcrop}" : ""
    MINLEN = ( params.trimmomatic.minlen > 0 ) ? "MINLEN:${params.trimmomatic.minlen}" : ""
    trimsteps = ILLUMINACLIP + ' ' + SLIDINGWINDOW + ' ' + LEADING + ' ' + TRAILING + ' ' + CROP + ' ' + HEADCROP + ' ' + MINLEN
    """
    trimmomatic PE -threads ${NumThreads} \
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
    input:
    set val(sampleName), file(readsFiles) from TrimmedReads

    output:
    tuple sampleName, file("${sampleName}.kraken"), file("${sampleName}.krpt") into KrakenFile
    tuple sampleName, file("${sampleName}.krpt") into KrakenVisuals

    script:
    quickflag = params.dev ? '--quick' : ''
    """
    kraken2 --db ${params.kraken.db} --threads ${NumThreads} --paired ${quickflag} \
        --report "${sampleName}.krpt" \
        --output "${sampleName}.kraken" \
        ${readsFiles}
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
    input:
    file(readsFiles) from PreKrakenReads
    set val(sampleName), file(krakenFile), file(krakenReport) from KrakenFile

    output:
    tuple sampleName, file("${sampleName}_filtered_{R1,R2}.fastq.gz") into ReadsForRay
    tuple sampleName, file("${sampleName}_filtered_{R1,R2}.fastq.gz") into ReadsForIVA
    tuple sampleName, file("${sampleName}_filtered_{R1,R2}.fastq.gz") into ReadsForMetaVelvet
    tuple sampleName, file("${sampleName}_filtered_{R1,R2}.fastq.gz") into ReadsForAbyss
    tuple sampleName, file("${sampleName}_filtered_{R1,R2}.fastq.gz") into ReadsForTrinity
    file("${sampleName}_filtered_{R1,R2}.fastq.gz") into CompressedReadsForRemapping

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

// Assemble into contigs using Ray
process ray {
    input:
    set val(sampleName), file(readsFiles) from ReadsForRay

    output:
    tuple val(sampleName), val(assembler), 'Contigs.fasta' into RayContigsForBlast
    tuple val(sampleName), val(assembler), file('Contigs.fasta'), file(readsFiles) into RayContigsForRemapping

    script:
    // Prevent command line from clobbering preset values
    kmerlength               = params.ray.kmerlength ?: 21
    minimumSeedLength        = params.ray.minimumSeedLength ?: 100
    minimumContigLength      = params.ray.minimumContigLength ?: 100
    maximumSeedCoverageDepth = params.ray.maximumSeedCoverageDepth ?: 4294967295
    minimumSeedCoverageDepth = params.ray.minimumSeedCoverageDepth ?: 0

    // Export the assembler for future combined steps
    assembler = 'ray'
    """
    mpiexec -n ${NumThreads} Ray -k ${kmerlength} \
        -minimumSeedLength ${minimumSeedLength} \
        -minimumContigLength ${minimumContigLength} \
        -maximumSeedCoverageDepth ${maximumSeedCoverageDepth} \
        -minimumSeedCoverageDepth ${minimumSeedCoverageDepth} -p ${readsFiles}
    mv RayOutput/Contigs.fasta .
    """
}

// Assemble into contigs using iva
process iva {
    input:
    set val(sampleName), file(readsFiles) from ReadsForIVA

    output:
    tuple val(sampleName), val(assembler), file('contigs.fasta') into IVAContigsForBlast
    tuple val(sampleName), val(assembler), file('contigs.fasta'), file(readsFiles) into IVAContigsForRemapping

    script:
    // Prevent parameter clobbering
    k                  = params.iva.k ?: 19
    s                  = params.iva.s ?: 11
    y                  = params.iva.y ?: 0.5
    ctg_first_trim     = params.iva.ctg_first_trim ?: 25
    ctg_iter_trim      = params.iva.ctg_iter_trim ?: 10
    ext_min_cov        = params.iva.ext_min_cov ?: 10
    ext_min_ratio      = params.iva.ext_min_ratio ?: 4
    ext_max_bases      = params.iva.ext_max_bases ?: 100
    ext_min_clip       = params.iva.ext_min_clip ?: 3
    max_contigs        = params.iva.max_contigs ?: 50
    seed_min_kmer_cov  = params.iva.seed_min_kmer_cov ?: 25
    seed_max_kmer_cov  = params.iva.seed_max_kmer_cov ?: 1000000
    seed_ext_max_bases = params.iva.seed_ext_max_bases ?: 50
    seed_ext_min_cov   = params.iva.seed_ext_min_cov ?: 10
    seed_ext_min_ratio = params.iva.seed_ext_min_ratio ?: 4
    max_insert         = params.iva.max_insert ?: 800
    assembler = 'iva'
    """
    iva -t ${NumThreads} -k ${k} -s ${s} -y ${y} \
        --ctg_first_trim ${ctg_first_trim} \
        --ctg_iter_trim ${ctg_iter_trim} \
        --ext_min_cov ${ext_min_cov} \
        --ext_min_ratio ${ext_min_ratio} \
        --ext_max_bases ${ext_max_bases} \
        --ext_min_clip ${ext_min_clip} \
        --max_contigs ${max_contigs} \
        --seed_min_kmer_cov ${seed_min_kmer_cov} \
        --seed_max_kmer_cov ${seed_max_kmer_cov} \
        --seed_ext_max_bases ${seed_ext_max_bases} \
        --seed_ext_min_cov ${seed_ext_min_cov} \
        --seed_ext_min_ratio ${seed_ext_min_ratio} \
        --max_insert ${max_insert} \
        -f ${readsFiles[0]} -r ${readsFiles[1]} out
    mv out/contigs.fasta .
    """
}

// Assemble using MetaVelvet
process metavelvet {
    input:
    set val(sampleName), file(readsFiles) from ReadsForMetaVelvet

    output:
    tuple val(sampleName), val(assembler), 'contigs.fa' into MetaVelvetContigsForBlast
    tuple val(sampleName), val(assembler), file('contigs.fa'), file(readsFiles) into MetaVelvetContigsForRemapping

    script:
    assembler = 'metavelvet'
    """
    velveth out 51 -fastq.gz -longPaired -separate ${readsFiles}
    velvetg out -exp_cov auto -ins_length 260
    meta-velvetg out
    mv out/contigs.fa .
    """
}

// Assemble using Abyss
process abyss {
    input:
    set val(sampleName), file(readsFiles) from ReadsForAbyss

    output:
    tuple val(sampleName), val(assembler), 'contigs.fa' into AbyssContigsForBlast
    tuple val(sampleName), val(assembler), file('contigs.fa'), file(readsFiles) into AbyssContigsForRemapping

    script:
    assembler = 'abyss'
    """
    abyss-pe name=${sampleName} k=21 in="${readsFiles}"
    cp ${sampleName}-contigs.fa contigs.fa
    """
}

// Assemble using Trinity
process trinity {
    input:
    set val(sampleName), file(readsFiles) from ReadsForTrinity

    output:
    tuple val(sampleName), val(assembler), 'Trinity.fasta' into TrinityContigsForBlast
    tuple val(sampleName), val(assembler), file('Trinity.fasta'), file(readsFiles) into TrinityContigsForRemapping

    script:
    assembler = 'abyss'
    """
    Trinity --seqType fq --left ${readsFiles[0]} --right ${readsFiles[1]} --CPU ${NumThreads} --max_memory 10G --output out
    mv out/Trinity.fasta .
    """
}

// Remap contigs using BWA
process bwa {
    input:
    set val(sampleName), val(assembler), file(contigs), file(readsFiles) from RayContigsForRemapping.concat(IVAContigsForRemapping, MetaVelvetContigsForRemapping, AbyssContigsForRemapping, TrinityContigsForRemapping)

    output:
    tuple val(sampleName), val(assembler), file(contigs), file("${sampleName}_${assembler}.sam") into RemappedReads

    script:
    // Fix clobbered parameters
    n = params.bwa.n ?: 0.04
    o = params.bwa.o ?: 1
    e = params.bwa.e ?: -1
    d = params.bwa.d ?: 16
    i = params.bwa.i ?: 5
    k = params.bwa.k ?: 2
    M = params.bwa.M ?: 3
    O = params.bwa.O ?: 4
    E = params.bwa.E ?: 4
    """
    cp ${readsFiles[0]} read1.fastq.gz
    cp ${readsFiles[1]} read2.fastq.gz
    gunzip read1.fastq.gz read2.fastq.gz
    bwa index ${contigs}
    bwa aln \
        -n ${n} \
        -o ${o} \
        -e ${e} \
        -d ${d} \
        -i ${i} \
        -k ${k} \
        -M ${M} \
        -O ${O} \
        -E ${E} \
        -t ${NumThreads} ${contigs} read1.fastq > ${sampleName}.1.sai
    bwa aln \
        -n ${n} \
        -o ${o} \
        -e ${e} \
        -d ${d} \
        -i ${i} \
        -k ${k} \
        -M ${M} \
        -O ${O} \
        -E ${E} \
        -t ${NumThreads} ${contigs} read2.fastq > ${sampleName}.2.sai
    bwa sampe ${contigs} \
        ${sampleName}.1.sai ${sampleName}.2.sai \
        ${readsFiles} > ${sampleName}_${assembler}.sam
    """
}

// Sort and compress the sam files for visualization
process sortsam {
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
    // This is a final step: publish it
    publishDir OutFolder, mode: 'copy'

    // Blast needs to happen on all contigs from all assemblers, and both
    // blastn and blastx needs to be applied to all contigs
    input:
    set val(sampleName), val(assembler), file(readsFiles) from RayContigsForBlast.concat(IVAContigsForBlast, MetaVelvetContigsForBlast, AbyssContigsForBlast, TrinityContigsForBlast)
    each program from BlastAlgorithms

    output:
    file("${RunName}_${sampleName}_${assembler}.${program}.tsv")

    script:
    // Prevent parameter clobbering
    max_hsps       = params.blast.max_hsps       ?: 10
    num_alignments = params.blast.num_alignments ?: 5
    outfmt         = params.blast.outfmt         ?: '"6 qseqid stitle sgi staxid ssciname scomname score bitscore qcovs evalue pident length slen saccver mismatch gapopen qstart qend sstart send"'
    evalue         = params.blast.evalue         ?: 1e-5

    // Pick the faster algorithm if this is a development cycle, otherwise
    // the titular program is also the name of the algorithm
    algorithm = program
    if ( params.dev ) {
        algorithm = ( program == 'blastn' ) ? 'megablast' : 'blastx-fast'
        max_hsps = 1
        num_alignments = 1
        evalue = 1e-50
    }

    // Switch which database to read from
    dbExtension = (program == 'blastn') ? 'nt' : 'nr'

    // Squash the filename into a single variable
    outFile = "${RunName}_${sampleName}_${assembler}.${program}.tsv"
    """
    echo "Sequence ID\tDescription\tGI\tTaxonomy ID\tScientific Name\tCommon Name\tRaw score\tBit score\tQuery Coverage\tE value\tPercent identical\tSubject length\tAlignment length\tAccession\tMismatches\tGap openings\tStart of alignment in query\tEnd of alignment in query\tStart of alignment in subject\tEnd of alignment in subject" > ${outFile}
    ${program} -query ${readsFiles} \
        -db ${params.blast.db}/${dbExtension} \
        -max_hsps ${max_hsps} \
        -num_alignments ${num_alignments} \
        -outfmt ${outfmt} \
        -evalue ${evalue} \
        -num_threads ${NumThreads} \
        -task ${algorithm} >> ${outFile}
    """
}
