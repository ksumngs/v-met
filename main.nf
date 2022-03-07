#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { KRAKEN2 } from './modules/ksumngs/nf-modules/kraken2/main.nf'
include { KRAKEN2_DBPREPARATION } from './modules/local/kraken2/dbpreparation.nf'
include { KRAKENTOOLS_EXTRACT } from './modules/ksumngs/nf-modules/krakentools/extract/main.nf'
include { QC } from './subworkflows/qc.nf'
include { READS_INGEST } from './subworkflows/ingest.nf'
include { TRIMMING } from './subworkflows/trimming.nf'
include { cowsay } from './lib/cowsay.nf'
include { vmet_logo } from './lib/logo.nf'

cowsay(vmet_logo())

if (params.help) {
    log.info(
        """\
        v-met - A bare-bones, ridiculously simple metagenomics pipeline for viruses

        Usage:

            nextflow run ksumngs/v-met

        Options:

            --input             Relative or absolute path to directory containing
                                gzipped fastq files or a TSV samplesheet
                                    type: path, default: .

            --platform          Type of reads to process. Options are 'illumina' and
                                'nanopore'
                                    type: string, default: none

            --kraken2_db        Kraken2-compatible database for classifying reads
                                    type: path, default: none

            --blast_db          Folder containing an NCBI BLAST NT database
                                    type: path, default: none

            --blast_target      Which reads to BLAST. Possible options are
                                'none'          Disable BLAST
                                'all'           BLAST all reads
                                'unclassified'  BLAST only Kraken's unclassified reads
                                'classified'    BLAST only Kraken's classified reads
                                - or -
                                A space-separated list of NCBI taxids to keep reads
                                that were classified as matching those taxids
                                Defaults to keeping unclassified and viral reads.
                                    type: string, default: '0 10239'

        For more information on usage and parameters, visit the website at
            https://ksumngs.github.io/v-met
        """.stripIndent()
        )

    exit 0
}

// Verify that the sequencing platform is specified
if (params.platform != 'illumina' && params.platform != 'nanopore') {
    log.error "ERROR: --platform <illumina,nanopore> must be specified"
    exit 1
}

// Verify that a kraken database was specified. We'll verify that it exists later
if (!params.kraken2_db) {
    log.error "ERROR: --kraken2_db must be specified"
    exit 1
}

// Verify that a BLAST database is specified if it is required
if ((!params.skip_blast && params.blast_target != 'none') && !params.blast_db) {
    log.error "ERROR: --blast_db must be specified, or else --blast_target 'none' must be passed"
    exit 1
}

log.info(
    """\
    Input folder:           ${params.input}
    Sequencing platform:    ${params.platform}
    Kraken2 Database:       ${params.kraken2_db}
    BLAST Database:         ${params.blast_db}
    Reads to BLAST:         '${params.blast_target}'
    Output folder           ${params.outdir}
    Diagnostics folder:     ${params.tracedir}
    """.stripIndent()
)

workflow {
    LogFiles = Channel.empty()
    VersionFiles = Channel.empty()

    // Bring in the reads files
    READS_INGEST()
    RawReads = READS_INGEST.out.sample_info
    VersionFiles = VersionFiles.mix(READS_INGEST.out.versions)

    // Check out the read quality
    if (!params.skip_qc) {
        QC(RawReads)
        LogFiles = LogFiles.mix(QC.out.report)
        VersionFiles = VersionFiles.mix(QC.out.versions)
    }

    // Trim the reads
    if (!params.skip_trimming) {
        TRIMMING(RawReads)
        TRIMMING.out.fastq.set{ TrimmedReads }
        LogFiles = LogFiles.mix(TRIMMING.out.log_out)
        VersionFiles = VersionFiles.mix(TRIMMING.out.versions)
    }
    else {
        RawReads.set{ TrimmedReads }
    }

    // Create a file object out of the Kraken2 database
    KrakenDb = file("${params.kraken2_db}", checkIfExists: true)
    if (!KrakenDb.isDirectory()) {
        // The Kraken database is not a local directory, so we'll assume it is a
        // tarballed database (local or remote)
        KRAKEN2_DBPREPARATION(KrakenDb)
        KrakenDb = KRAKEN2_DBPREPARATION.out.db
        VersionFiles = VersionFiles.mix(KRAKEN2_DBPREPARATION.out.versions)
    }

    // Do Kraken on the reads
    KRAKEN2(TrimmedReads, KrakenDb)
    VersionFiles = VersionFiles.mix(KRAKEN2.out.versions)
    LogFiles = LogFiles.mix(KRAKEN2.out.kreport)

    KrakenReports = KRAKEN2.out.reports
    UnclassifiedReads = KRAKEN2.out.unclassified
    ClassifiedReads = KRAKEN2.out.classified

    kraken2krona(KrakenReports)
    KrakenKronas = kraken2krona.out

    krona(KrakenKronas.collect())

    if (!params.skip_blast && params.blast_target != 'none') {
        if (params.blast_target == 'all') {
            TrimmedReads | convert2fasta | blast
        }
        else if (params.blast_target == 'classified') {
            ClassifiedReads | convert2fasta | blast
        }
        else if (params.blast_target == 'unclassified') {
            UnclassifiedReads | convert2fasta | blast
        }
        else {
            KrakenReads = TrimmedReads.join(KrakenReports)
            KrakenReads | filterreads | convert2fasta | blast
        }
    }
}

process sample_rename {
    label 'process_low'

    input:
    tuple val(givenName), file(readsFiles)

    output:
    tuple val(sampleName), file("out/*.fastq.gz")

    script:
    sampleName = givenName.split('_')[0]
    if (params.ont) {
        """
        mkdir out
        mv ${readsFiles[0]} out/${sampleName}.fastq.gz
        """
    }
    else {
        """
        mkdir out
        mv ${readsFiles[0]} out/${sampleName}_R1.fastq.gz
        mv ${readsFiles[1]} out/${sampleName}_R2.fastq.gz
        """
    }
}

process nanofilt {
    label 'nanofilt'
    label 'process_low'

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    tuple val(sampleName), path("${sampleName}_trimmed.fastq.gz")

    script:
    minlenflag = ( params.trim_minlen > 0 )   ? "--length ${params.trim_minlen}"     : ''
    maxlenflag = ( params.trim_maxlen > 0 )   ? "--maxlength ${params.trim_maxlen}"  : ''
    qualflag   = ( params.trim_meanqual > 0 ) ? "--quality ${params.trim_meanqual}"  : ''
    mingcflag  = ( params.trim_mingc > 0 )    ? "--minGC ${params.trim_mingc}"       : ''
    maxgcflag  = ( params.trim_maxgc > 0 )    ? "--maxGC ${params.trim_maxgc}"       : ''
    headflag   = ( params.trim_headcrop > 0 ) ? "--headcrop ${params.trim_headcrop}" : ''
    tailflag   = ( params.trim_tailcrop > 0 ) ? "--tailcrop ${params.trim_tailcrop}" : ''
    optionflags = [
        minlenflag,
        maxlenflag,
        qualflag,
        mingcflag,
        maxgcflag,
        headflag,
        tailflag
    ].join(' ')
    """
    gunzip < ${readsFiles} | \
        NanoFilt --logfile ${sampleName}.nanofilt.log ${optionflags} | \
        gzip -9 > ${sampleName}_trimmed.fastq.gz
    """
}


// First trim, using Trimmomatic
// Trim Illumina reads
process trimmomatic {
    label 'trimmomatic'
    label 'process_medium'

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    tuple val(sampleName), path("*.fastq.gz")

    script:
    // Put together the trimmomatic parameters
    clipflag =
        ( !(params.trim_adapters.getClass() == Boolean || params.trim_adapters.allWhitespace) &&
        params.trim_mismatches > 0 && params.trim_pclip > 0 && params.trim_clip) ?
        "ILLUMINACLIP:/Trimmomatic-0.39/adapters/${params.trim_adapters}:${params.trim_mismatches}:${params.trim_pclip}:${params.trim_clip}" : ''
    winflag =
        ( params.trim_winsize > 0 && params.trim_winqual > 0 ) ? "SLIDINGWINDOW:${params.trim_winsize}:${params.trim_winqual}" : ""
    leadflag =
        ( params.trim_leading > 0 ) ? "LEADING:${params.trim_leading}" : ""
    trailflag =
        ( params.trim_trailing > 0 ) ? "TRAILING:${params.trim_trailing}" : ""
    cropflag =
        ( params.trim_tailcrop > 0 ) ? "CROP:${params.trim_crop}" : ""
    headflag =
        ( params.trim_headcrop > 0 ) ? "HEADCROP:${params.trim_headcrop}" : ""
    minlenflag =
        ( params.trim_minlen > 0 ) ? "MINLEN:${params.trim_minlen}" : ""
    trimsteps = [
        clipflag,
        winflag,
        leadflag,
        trailflag,
        cropflag,
        headflag,
        minlenflag
    ].join(' ')
    """
    trimmomatic PE -threads ${task.cpus} \
        ${readsFiles} \
        ${sampleName}_trimmed_R1.fastq.gz \
        /dev/null \
        ${sampleName}_trimmed_R2.fastq.gz \
        /dev/null \
        ${trimsteps} 2> ${sampleName}.trimmomatic.log
    """
}

// Classify reads using Kraken
process kraken {
    label 'kraken'
    label 'process_high_memory'

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    tuple val(sampleName), file("${sampleName}.kraken"), file("${sampleName}.kreport"), emit: reports
    tuple val(sampleName), path("${sampleName}_classified*.fastq.gz"), emit: classified
    tuple val(sampleName), path("${sampleName}_unclassified*.fastq.gz"), emit: unclassified

    script:
    pairedflag       = params.pe ? '--paired' : ''
    classifiedflag   = params.pe ? "${sampleName}_classified#.fastq"   : "${sampleName}_classified.fastq"
    unclassifiedflag = params.pe ? "${sampleName}_unclassified#.fastq" : "${sampleName}_unclassified.fastq"
    """
    kraken2 --db ${params.kraken2_db} --threads ${task.cpus} ${pairedflag} \
        --report "${sampleName}.kreport" \
        --output "${sampleName}.kraken" \
        --classified-out ${classifiedflag} \
        --unclassified-out ${unclassifiedflag} \
        ${readsFiles}
    gzip -9 *.fastq
    """
}

// Convert the kraken reports to krona input files using KrakenTools then
// prettify them using unix tools
process kraken2krona {
    label 'krakentools'
    label 'process_low'

    input:
    tuple val(sampleName), file(krakenFile), file(krakenReport)

    output:
    file "${sampleName}.krona"

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
    label 'krona'
    label 'process_low'

    // This is a final step: publish it
    publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

    input:
    file '*'

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
    label 'krakentools'
    label 'process_low'

    input:
    tuple val(sampleName), file(readsFiles), file(krakenFile), file(krakenReport)

    output:
    tuple val(sampleName), file("${sampleName}_filtered_{R1,R2}.fastq.gz")

    // Although I haven't seen it documented anywhere, 0 is unclassified reads
    // and 10239 is viral reads
    script:
    """
    extract_kraken_reads.py -k ${krakenFile} \
        -s1 ${readsFiles[0]} -s2 ${readsFiles[1]} \
        -r ${krakenReport} \
        -t ${params.blast_target} --include-children \
        --fastq-output \
        -o ${sampleName}_filtered_R1.fastq -o2 ${sampleName}_filtered_R2.fastq
    gzip ${sampleName}_filtered_{R1,R2}.fastq
    """

}

// Convert the viral and unclassified reads to fasta for blast
process convert2fasta {
    label 'seqtk'
    label 'process_low'

    input:
    tuple val(sampleName), file(readsFiles)

    output:
    tuple val(sampleName), file("${sampleName}.fasta")

    script:
    """
    seqtk mergepe ${readsFiles} | seqtk seq -a - > ${sampleName}.fasta
    """
}

// Blast contigs
process blast {
    label 'blast'
    label 'process_medium'

    // This is a final step: publish it
    publishDir "${params.outdir}", mode: "${params.publish_dir_mode}"

    // Blast needs to happen on all contigs from all assemblers, and both
    // blastn and blastx needs to be applied to all contigs
    input:
    tuple val(sampleName), file(readsFiles)

    output:
    file("${sampleName}.blast.tsv")

    script:
    // Separate parameters from script
    max_hsps       = 10
    num_alignments = 5
    outfmt         = '"6 qseqid stitle sgi staxid ssciname scomname score bitscore qcovs evalue pident length slen saccver mismatch gapopen qstart qend sstart send"'
    evalue         = 1e-5

    // Squash the filename into a single variable
    outFile = "${sampleName}.blast.tsv"
    """
    echo "Sequence ID\tDescription\tGI\tTaxonomy ID\tScientific Name\tCommon Name\tRaw score\tBit score\tQuery Coverage\tE value\tPercent identical\tSubject length\tAlignment length\tAccession\tMismatches\tGap openings\tStart of alignment in query\tEnd of alignment in query\tStart of alignment in subject\tEnd of alignment in subject" > ${outFile}
    blastn -query ${readsFiles} \
        -db ${params.blast_db}/nt \
        -max_hsps ${max_hsps} \
        -num_alignments ${num_alignments} \
        -outfmt ${outfmt} \
        -evalue ${evalue} \
        -num_threads ${task.cpus} >> ${outFile}
    """
}
