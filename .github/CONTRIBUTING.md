# v-met: Contributing Guidelines

Hi there! Many thanks for taking an interest in improving v-met.

We try to manage the required tasks for v-met using GitHub issues, you probably
came to this page when creating one. Please use the pre-filled template to save
time.

However, don't be put off by this template - other more general issues and
suggestions are welcome! Contributions to the code are even more welcome ;)

> v-met comforms to pipeline standards made by the [nf-core
> group](https://nf-co.re). These guidelines are a mix of general GitHub
> workflow, nf-core compliance, and our own way of doing things.

## Contribution workflow

If you'd like to write some code for v-met, the standard workflow is as
follows:

1. Check that there isn't already an issue about your idea in the
   [ksumngs/v-met issues](https://github.com/ksumngs/v-met/issues) to avoid
   duplicating work
    - If there isn't one already, please create one so that others know you're
      working on this
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo)
   the [ksumngs/v-met repository](https://github.com/ksumngs/v-met) to
   your GitHub account
3. Make the necessary changes / additions within your forked repository
   following [Pipeline conventions](#pipeline-contribution-conventions)
4. Use `nf-core schema build` and add any new parameters to the pipeline JSON
   schema (requires [nf-core tools](https://github.com/nf-core/tools) >= 1.10).
5. Submit a Pull Request against the `dev` branch and wait for the code to be
   reviewed and merged

If you're not used to this workflow with git, you can start with some [docs from
GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests)
or even their [excellent `git` resources](https://try.github.io/).

## Tests

When you create a pull request with changes, [GitHub
Actions](https://github.com/features/actions) will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing,
though of course we can help out before then.

There are typically two types of tests that run:

### Lint tests

`nf-core` has a [set of guidelines](https://nf-co.re/developers/guidelines)
which we strive to adhere to. To enforce these and ensure that all pipelines
stay in sync, there is a helper tool which runs checks on the pipeline code.
This is in the [nf-core/tools repository](https://github.com/nf-core/tools) and
once installed can be run locally with the `nf-core lint <pipeline-directory>`
command.

> :bulb: To install nf-core tools on Linux or Mac, run the following commands
>
> ```bash
> python3 -m pip install --user pipx
> python3 -m pipx ensurepath
> pipx install nf-core
> ```

If any failures or warnings are encountered, please follow the listed URL for
more documentation.

### Pipeline tests

We use the test dataset over at
<https://github.com/ksumngs/nf-test-datasets/tree/v-met> for testing. `GitHub
Actions` then runs the pipeline on this data to ensure that it exits
successfully. If there are any failures then the automated tests fail. These
tests are run both with the latest available version of `Nextflow` and also the
minimum required version that is stated in the pipeline code.

## Patch

:warning: Only in the unlikely and regretful event of a release happening with a
bug.

- On your own fork, make a new branch `patch` based on `upstream/master`.
- Fix the bug, and bump version (X.Y.Z+1).
- A PR should be made on `master` from patch to directly this particular bug.

## Getting help

For further information/help, please consult the [v-met
documentation](https://ksumngs.github.io/v-met).

## Pipeline contribution conventions

To make the v-met code and processing logic more understandable for new
contributors and to ensure quality, we semi-standardize the way the code and
other contributions are written.

### Adding a new step

If you wish to contribute a new step, please use the following coding standards:

1. Define the corresponding input channel into your new process from the
   expected previous process channel
2. Write the process block (see below).
3. Define the output channel if needed (see below).
4. Add any new parameters to `nextflow.config` with a default (see below).
5. Add any new parameters to `nextflow_schema.json` with help text (via the
   `nf-core schema build` tool).
6. Add sanity checks and validation for all relevant parameters.
7. Perform local tests to validate that the new code works as expected.
8. If applicable, add a new test command in `.github/workflow/ci.yml`.
9. Update MultiQC config `assets/multiqc_config.yaml` so relevant suffixes, file
   name clean up and module plots are in the appropriate order. If applicable,
   add a [MultiQC](https://https://multiqc.info/) module.
10. Add a description of the output files and if relevant any appropriate images
    from the MultiQC report to `docs/output.md`.

### Default values

Parameters should be initialised / defined with default values in
`nextflow.config` under the `params` scope.

Once there, use `nf-core schema build` to add to `nextflow_schema.json`.

### Default processes resource requirements

Every process definition should be specified with `withLabel:` selectors so they
will request the proper amount of resources. A nf-core standard set of labels
that should be followed where possible can be seen in the [nf-core pipeline
template](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/conf/base.config),
which has the default process as a single core-process, and then different
levels of multi-core configurations for increasingly large memory requirements
defined with standardised labels.

The process resources can be passed on to the tool dynamically within the
process with the `${task.cpu}` and `${task.memory}` variables in the `script:`
block.

### Naming schemes

There are only two hard things in computer science: naming and cache
invalidation. With regards to the later, in case of failure remove the `-resume`
flag. With regards to the former, we name things slightly different from
nf-core's guidelines. Specifically, we use these [naming
guidelines](https://gist.github.com/MillironX/bd9606623b3ccfdfb72d77e2bd3dc213#naming):

> - Use `UPPER_CASE_WITH_UNDERSCORE_SEPARATORS` for Process and Workflow names
> - Use `CapitalizeEveryWordWithoutSeparators` for Channel and Global variable
>   names
>   - Exception: output channels
> - Use `snake_case_underscore_separators` for parameter and process label
>   names, and input and output variable names


### Nextflow version bumping

If you are using a new feature from core Nextflow, you may bump the minimum
required version of nextflow in the pipeline with: `nf-core bump-version
--nextflow . [min-nf-version]`

### Images and figures

For overview images and other documents we follow the nf-core [style guidelines
and examples](https://nf-co.re/developers/design_guidelines).
