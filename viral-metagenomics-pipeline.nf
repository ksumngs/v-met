#!/usr/bin/env nextflow

params.R1
params.R2
params.threads
params.out

Channel
    .fromPath(params.R1)
    .set{ ForwardRead }
Channel
    .fromPath(params.R2)
    .set{ ReverseRead }

NumThreads = params.threads

process trim {
    input:
        path ForwardRead
        path ReverseRead
        val NumThreads

    output:
        file 'paired_out_R1.fastq.gz' into Paired1
        file 'unpaired_out_R1.fastq.gz' into Unpaired1
        file 'paired_out_R2.fastq.gz' into Paired2
        file 'unpaired_out_R2.fastq.gz' into Unpaired2

        shell:
        '''
            trimmomatic PE -threads !{NumThreads} \
                !{ForwardRead} \
                !{ReverseRead} \
                paired_out_R1.fastq.gz \
                unpaired_out_R1.fastq.gz \
                paired_out_R2.fastq.gz \
                unpaired_out_R2.fastq.gz \
                ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        '''
}

Paired1
    .collectFile()
Paired2
    .collectFile()
Unpaired1
    .collectFile()
Unpaired2
    .collectFile()
