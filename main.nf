#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BLAST_ADDHEADER } from './modules/local/blast/addheader.nf'
include { BLAST_BLASTN } from './modules/nf-core/modules/blast/blastn/main.nf'
include { BLAST_DBPREPARATION } from './modules/local/blast/dbpreparation.nf'
include { CAT_FASTQ as UNZIP_FASTA } from './modules/ksumngs/nf-modules/cat/fastq/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/modules/custom/dumpsoftwareversions/main.nf'
include { KRAKEN2 } from './modules/ksumngs/nf-modules/kraken2/main.nf'
include { KRAKEN2_DBPREPARATION } from './modules/local/kraken2/dbpreparation.nf'
include { KRAKENTOOLS_EXTRACT } from './modules/ksumngs/nf-modules/krakentools/extract/main.nf'
include { KRAKENTOOLS_KREPORT2KRONA } from './modules/ksumngs/nf-modules/krakentools/kreport2krona/main.nf'
include { KRONA_IMPORTTEXT } from './modules/ksumngs/nf-modules/krona/importtext/main.nf'
include { MULTIQC } from './modules/nf-core/modules/multiqc/main.nf'
include { QC } from './subworkflows/qc.nf'
include { READS_INGEST } from './subworkflows/ingest.nf'
include { SEQTK_SEQ } from './modules/nf-core/modules/seqtk/seq/main.nf'
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

    // Make a Krona chart
    KRAKENTOOLS_KREPORT2KRONA(KRAKEN2.out.kreport)
    KRONA_IMPORTTEXT(
        KRAKENTOOLS_KREPORT2KRONA.out.krona
            .map{ it.drop(1) }
            .collect()
        )
    VersionFiles = VersionFiles.mix(KRAKENTOOLS_KREPORT2KRONA.out.versions)
    VersionFiles = VersionFiles.mix(KRONA_IMPORTTEXT.out.versions)

    // BLAST reads, if requested
    if (!params.skip_blast && params.blast_target != 'none') {
        // Remove the BLAST target reads into their own channel
        if (params.blast_target == 'all') {
            TrimmedReads.set{ FilteredReads }
        }
        else if (params.blast_target == 'classified') {
            KRAKEN2.out.classified.set{ FilteredReads }
        }
        else if (params.blast_target == 'unclassified') {
            KRAKEN2.out.unclassified.set{ FilteredReads }
        }
        else {
            TrimmedReads
                .join(KRAKEN2.out.kraken)
                .join(KRAKEN2.out.kreport)
                .set{ ReadPlusReports }
            KRAKENTOOLS_EXTRACT(ReadPlusReports, "${params.blast_target}")
            KRAKENTOOLS_EXTRACT.out.fastq.set{ FilteredReads }
            VersionFiles = VersionFiles.mix(KRAKENTOOLS_EXTRACT.out.versions)
        }

        // Convert the reads to FASTA to be a suitable BLAST query
        SEQTK_SEQ(FilteredReads)
        VersionFiles = VersionFiles.mix(SEQTK_SEQ.out.versions)
        UNZIP_FASTA(
            SEQTK_SEQ.out.fastx
                .map{[
                    ['id': it[0].id, 'single_end': true],
                    it[1]
                ]},
            false
        )
        VersionFiles = VersionFiles.mix(UNZIP_FASTA.out.versions)

        // Create a file object out of the BLAST database
        BlastDb = file("${params.blast_db}", checkIfExists: true)
        if (!BlastDb.isDirectory()) {
            // The BLAST database is not a local directory, so we'll assume it is a
            // tarballed database (local or remote)
            BLAST_DBPREPARATION(BlastDb)
            BlastDb = BLAST_DBPREPARATION.out.db
            VersionFiles = VersionFiles.mix(BLAST_DBPREPARATION.out.versions)
        }

        // BLAST reads
        BLAST_BLASTN(UNZIP_FASTA.out.reads, BlastDb)
        VersionFiles = VersionFiles.mix(BLAST_BLASTN.out.versions)

        // Add header to BLAST output
        BLAST_ADDHEADER(
            BLAST_BLASTN.out.txt,
            file("${workflow.projectDir}/assets/headers/blast_outfmt6_header.txt", type: 'file', checkIfExists: true)
        )
    }

    // Collect the logs
    LogFiles = LogFiles
        .map{ (it instanceof Path) ? it : it.drop(1) }
        .mix(Channel.of(file("${workflow.projectDir}/assets/multiqc_config.yaml")))
        .flatten()
        .collect()

    // Create a MultiQC report
    MULTIQC(LogFiles)
    VersionFiles = VersionFiles.mix(MULTIQC.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS(
        VersionFiles
            .unique()
            .collectFile(name: 'collated_versions.yml')
    )
}
