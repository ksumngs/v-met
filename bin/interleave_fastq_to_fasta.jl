#!/usr/bin/env julia

# Install the FASTX package
using Pkg
required_packages = ["FASTX", "CodecZlib"]
for pkg in required_packages
    if !(pkg ∈ keys(Pkg.installed()))
        Pkg.add(pkg)
    end
end

# Import packages
using FASTX
using CodecZlib

# Import arguments
fwdread = ARGS[1]
revread = ARGS[2]
combinedreads = ARGS[3]

# Prepare the IO streams
fwdreader = FASTQ.Reader(GzipDecompressorStream(open(fwdread, "r")))
revreader = FASTQ.Reader(GzipDecompressorStream(open(revread, "r")))
cmbwriter = FASTA.Writer(open(combinedreads, "w"))

# Copy the records in, alternating between forward and reverse reads
while !eof(fwdreader)
    fastqfwdrecord = read(fwdreader)
    fastafwdrecord = FASTA.Record(
        identifier(fastqfwdrecord),
        description(fastqfwdrecord),
        sequence(fastqfwdrecord)
    )
    fastqrevrecord = read(revreader)
    fastarevrecord = FASTA.Record(
        identifier(fastqrevrecord),
        description(fastqrevrecord),
        sequence(fastqrevrecord)
    )
    write(cmbwriter, fastafwdrecord)
    write(cmbwriter, fastarevrecord)
end

# Close the streams
close(fwdreader)
close(revreader)
close(cmbwriter)
