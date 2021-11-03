Usage
=====

Basic jist::

    nextflow run ksumngs/v-met -profile <docker,podman,singularity> --platform <pe,ont> --kraken2_db <kraken2_db> --blast_db <blast_db>

Input Preparation
-----------------

Input is made up of raw NGS reads in gzipped fastq format. There must be only
one file per sample for Nanopore reads, and one pair of files (two files) for
Illumina paired-end reads. Paired-end files must have either the characters
``_1`` and ``_2`` or ``_R1`` and ``_R2`` in the filename to be idenfied as a
pair. Note that this pipeline does not support Illumina single-end reads. The
sample name is taken as the all the characters *before* the first underscore in
the file name.

For example, a folder with the following files

::

    .
    ├── pig-serum_S4_L001_R1_001.fastq.gz
    ├── pig-serum_S4_L001_R2_001.fastq.gz
    ├── pig-feces_S5_L001_R1_001.fastq.gz
    ├── pig-feces_S5_L001_R2_001.fastq.gz
    ├── mosquito1_S6_L001_R1_001.fastq.gz
    ├── mosquito1_S6_L001_R2_001.fastq.gz
    ├── mosquito2_S7_L001_R1_001.fastq.gz
    └── mosquito2_S7_L001_R2_001.fastq.gz

will produce a sample list of

* pig-serum
* pig-feces
* mosquito1
* mosquito2

This file naming system is very common for Illumina output, however the default
Oxford Nanopore naming scheme breaks this in several ways:

1. The fastq files are (often) not gzipped
2. The fastq files are stored in separate directories per sample
3. The fastq files start with the flowcell id as the beginning of the filename
4. There are often many fastq files per sample

For example, MinION output folders often are structured like the following

::

    .
    ├── fastq_fail
    │   ├── barcode04
    |   |   ├── FAP01234_fail_barcode04_abcd1234_0.fastq
    |   |   └── FAP01234_fail_barcode04_abcd1234_1.fastq
    |   ├── barcode05
    |   |   ├── FAP01234_fail_barcode05_abcd1234_0.fastq
    |   |   └── FAP01234_fail_barcode05_abcd1234_1.fastq
    │   ├── barcode06
    |   |   └── FAP01234_fail_barcode06_abcd1234_0.fastq
    |   └── barcode07
    |       ├── FAP01234_fail_barcode07_abcd1234_0.fastq
    |       └── FAP01234_fail_barcode07_abcd1234_1.fastq
    └── fastq_pass
        ├── barcode04
        |   ├── FAP01234_pass_barcode04_abcd1234_0.fastq
        |   └── FAP01234_pass_barcode04_abcd1234_1.fastq
        ├── barcode05
        |   ├── FAP01234_pass_barcode05_abcd1234_0.fastq
        |   └── FAP01234_pass_barcode05_abcd1234_1.fastq
        ├── barcode06
        |   └── FAP01234_pass_barcode06_abcd1234_0.fastq
        └── barcode07
            ├── FAP01234_pass_barcode07_abcd1234_0.fastq
            └── FAP01234_pass_barcode07_abcd1234_1.fastq


(Non-fastq files and directories excluded for brevity.)

One way of remedying these difficulties is to create a tab-separated file with
the sample ID and the barcode (with leading zeros), such as

========= =======
sample    barcode
========= =======
pig-serum 04
pig-feces 05
mosquito1 06
mosquito2 07
========= =======

Then the following bash script can restructure the files into the correct
format using GNU Parallel and pigz, and can handle a mix of gzipped and
uncompressed reads. Replace ``fastq*`` with ``fastq_pass`` in line 8 if you
wish to only consider those reads that passed Guppy's quality checks.

.. code-block:: bash

    #!/bin/bash
    SAMPLESHEET=/path/to/samplesheet
    MINION_DIR=/path/to/minion/output
    while read -r LINE; do
      echo "$LINE" | while IFS=$'\t' read -r SAMPLE BARCODE; do
        if [[ $SAMPLE != "sample" ]]; then
          mkdir $SAMPLE
          cp ${MINION_DIR}/fastq*/barcode${BARCODE}/*.fastq* $SAMPLE
          parallel gunzip ::: $SAMPLE/*.fastq.gz
          cat $SAMPLE/*.fastq > $SAMPLE.fastq
          pigz --best -p $(nproc) $SAMPLE.fastq
          rm -rf $SAMPLE
        fi
      done
    done < $SAMPLESHEET

Running this script on the example directory and sample sheet will yield the
following files

::

    .
    ├── pig-serum.fastq.gz
    ├── pig-feces.fastq.gz
    ├── mosquito1.fastq.gz
    └── mosquito2.fastq.gz

and will give the same sample names as given above for the Illumina files. Note
that fast5 files are not needed for Nanopore reads.

Profile Selection
-----------------

Profiles allow for unique combinations of settings within a Nextflow pipeline.
For the purposes of he v-met pipeline, they reconfigure the pipeline to run
on a particular container engine. Whichever engine you choose must be installed
and available (e.g. ``module load``) on each node that pipeline processes are
running on. The available options are

(none)
  Don't use a container engine. Requires that every tool and the right version
  of every tool be installed on each machine the pipeline is running on. Don't
  use this one.
docker
  Use `Docker <https://docker.com>`_ as the container engine. Note that Docker
  often requires root or nearly-root permissions that usually aren't available
  on HPCs, and has a weird license that might forbid its use in commercial
  settings. Works well on local machines, though.
podman
  Use `Podman <https://podman.io>`_ instead of Docker. Podman is similar enough
  to Docker that they can often be used interchangably, but doesn't require root
  permissions and has a free license. Some files might not be accessible via
  container on RHEL-based distros thanks to their particular SELinux
  implementation.
singularity
  **Recommended**

  Use the
  :doc:`Singularity <singularity:index>` container enginer.
  This engine was build with HPC use in mind, and doesn't require any special
  permissions to run. Singularity's downfall is how it will expose your home
  directory to the container, resulting in rare, but difficult to explain bugs
  when files conflict. Every effort has been made to minimize the effects of
  Singularity's file mounting in this pipeline.

To select a profile, you must pass the desired profile name to Nextflow's
``-profile`` flag. Note that this is a Nextflow flag, and not a pipeline flag,
so it is a single dash (``-profile``), not a double dash (``--profile``).

Mandatory Parameters
--------------------

See :doc:`the page on parameters <parameters>` for the complete lowdown on
parameters.

The pipeline is pretty much set up to run itself. So long as you have your input
reads formatted correctly, it doesn't need much input from you. These are the
bare minimum parameters that you must provide on the command-line for the
pipeline to complete. Note that there has been mixed success with placing these
parameters in a ``nextflow.config`` file, so keeping them on the command-line is
best.

--kraken2_db
  The path to a Kraken2 database. See :ref:`--kraken2_db`.
--blast_db
  The path to an NCBI nt BLAST database. See :ref:`--blast_db`.
--platform
  Must be set to ``ont`` or ``pe``, depending on the type of reads you are
  analyzing. See :ref:`--platform`.
-profile
  So, this isn't really a parameter, but the container engine needs to be set
  using Nextflow's ``-profile`` flag. See :ref:`Profile Selection`.

Setting up for HPC Job Schedulers
---------------------------------

v-met comes preconfigured for local use only. Yes, that's about as ridiculous as
it sounds. What's even more ridiculous is trying to make a configuration that
can easily be adapted to multiple HPCs and job-scheduler frameworks. There is a
compromise, however.

Process Labels
^^^^^^^^^^^^^^

Rather than provide hard-coded configurations that will certainly
break, there are several Nextflow 'labels,' that can be used for assigning
processes to specific node queues or partitions. The labels are

process_low
  For processes with low resource usage that take a short time
process_medium
  For processes with moderate resource usage that take a moderate amount of time
process_high
  For processes with high resource usage that take a moderately high amount of
  time
process_long
  For processes that take a long time
process_high_memory
  For processes that use a lot of memory
run_local
  For processes that have to be run on the login node. This label was created
  specially for processes that download resources from the internet on HPCs
  where compute nodes do not have internet access

Using a custom ``nextflow.config`` and these process labels, you can construct a
setup for your own HPC needs.

Example
^^^^^^^

As an example, here is a guide on how to set up a configuration for the
`USDA's SCINet Ceres cluster <https://scinet.usda.gov/guide/ceres/>`_, using the
publically available info on their website.

First, we see that Ceres uses SLURM and Singularity. Excellent.
Let's set Nextflow up to use SLURM:

.. code-block:: groovy

    process {
        executor = 'slurm'
    }

Some SLURM systems require an account to be passed with every job submission.
Let's add ours just to be safe.

.. code-block:: groovy

    process {
        executor       = 'slurm'
        clusterOptions = '--account=ksumngs'
    }

For this example, I don't think we'll need to do anything special with the low,
medium, and high processes, but let's make sure that the long and high memory
processes get submitted to partitions that can handle them.

.. code-block:: groovy

    process {
        executor       = 'slurm'
        clusterOptions = '--account=ksumngs'
        module         = 'singularity'
        withLabel: process_long {
            queue      = 'long'
        }
        withLabel: process_high_memory {
            queue      = 'mem'
        }
    }

Now, we can place this file in ``$HOME/.nextflow/nextflow.config``, and these
settings will be applied every time we run the pipeline.
