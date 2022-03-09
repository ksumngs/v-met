---
name: Bug Report
about: Something isn't working like it's supposed to
title: "[Bug]: "
labels: ["bug"]
---

> _Thanks for helping squash bugs in v-met!_ :smile:
> _This template is rather extensive._
> _Please fill out all that you can, but remember that this template_
> _is only a tool: if you feel that anything is not relevant, simply_
> _leave it blank._

## Behavior

### Expected behavior

> _Tell us what you expected to happen_

<!-- It should just work! -->

### Current behavior

> _Tell us what actually happened_

<!-- A bug happened! -->

## Steps to reproduce

### File(s) that cause bug

> _Please drag and drop (i.e. upload to the GitHub issue) the exact files that
> are causing the problem. If needed, you can zip/tarball them first._

<!--
- [sample.bam](#)
- [reference.fasta](#)
-->

### Command(s) used to run the pipeline

> _List the shell commands used to run the pipeline on the files you uploaded._

<!--
```bash
nextflow run ksumngs/v-met \
  -profile singularity      \
  --platform illumina       \
  --kraken2_db /databases/kraken2/nt
```
-->

## Pipeline Output

### Log file

> _Please drag and drop the `.nextflow.log` file from this pipeline run_

### Results folder

> _If results are generated (and the bug is related to invalid output), please
> zip/tarball the results folder and upload it_

## Additional Info

### Context

> _Tell us what you are trying to accomplish and/or how this bug affects users
> in the "real world."_

### Possible Solution

> _**(Highly optional, but highly useful.)**
> If you can figure out a design decision or a particular function call
> that is causing the bug, explain here.
> Even better, if you know a way to fix it, but don't want to open a pull
> request, then explain the fix here._

## Environment

- Hardware: <!-- [e.g. HPC, Desktop, Cloud...] -->
- Executor: <!-- [e.g. slurm, local, awsbatch...] -->
- OS: <!-- [e.g. CentOS Linux, macOS, Linux Mint...] -->
- Nextflow Version: <!-- [e.g. 21.04.0] -->
- Container Engine: <!-- [e.g. Docker, Singularity, Podman] -->
