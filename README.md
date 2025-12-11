<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-swgsrelate_logo_dark.png">
    <img alt="nf-core/swgsrelate" src="docs/images/nf-core-swgsrelate_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/swgsrelate)
[![GitHub Actions CI Status](https://github.com/nf-core/swgsrelate/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/swgsrelate/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/swgsrelate/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/swgsrelate/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/swgsrelate/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/swgsrelate)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23swgsrelate-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/swgsrelate)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/swgsrelate** is a bioinformatics pipeline that for estimating genetic relatedness from low-coverage whole-genome sequencing (sWGS) data. It performs read mapping, optional base quality score recalibration, variant calling with GATK and BCFtools, and downstream relatedness estimation using multiple complementary tools. For many non-model organisms, no high-confidence variant set is available. The pipeline provides an automated multi-round bootstrapping workflow to generate one. The resulting standardized outputs include genotype likelihood-based variant calls, filtered VCF files, and relatedness estimates from several independent algorithms, enabling robust inference even from very sparse sequencing data.

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

The pipeline can perform the following major processing stages:

1. Input parsing & metadata setup: Reads a CSV samplesheet describing the input FASTQ, SPRING or CRAM files and their read-group information.

2. Reference genome preparation:
   If not provieded, this step automatically generates:

- BWA-MEM2 index files
- FASTA index (.fai)
- Sequence dictionary (.dict)

3. Read alignment

- Aligns raw FASTQ reads to the reference genome.
- Produces sorted, indexed CRAM files with proper read-group annotations.

4. If a known variant set is provided, runs Base Quality Score Recalibration.

5. If no known variant set is available, the pipeline can generate one automatically via bootstrapping.

- Iteratively refines and stabilises the set of high-confidence SNPs for downstream use.

6. Runs Base Quality Score Recalibration if a known variant set was provided.

7. Variant calling (GATK HaplotypeCaller & BCFtools)

8. Performs joint variant discovery for all samples.

9. Combines GATK and BCFtools results using bcftools isec to produce a conservative, high-confidence set of variants.

10. Variant filtering & thinning and optional exclusion of specific scaffolds.

11. Relatedness estimation (multi-tool)
    Uses multiple complementary tools to increase robustness, depending on configuration:

- READv2 (ML-based relatedness estimation for low-coverage data)
- BREADR (R-based Bayesian relatedness inference)
- NGSrelate/ANGSD (likelihood-based estimation directly from genotype likelihoods)

12. MultiQC reporting: Aggregates quality metrics across all workflow stages into a single interactive report.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):-->

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz

```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end). Alternatively, the samplesheet can be filled with fastq files encoded in SPRING format, or the preprocessing steps can be skipped entirely when BAM or CRAM files are provided.

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run nf-core/swgsrelate \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --fasta <REFGENOME>
   --bootstrapping_rounds 1\
   --outdir <OUTDIR>
```

> **Note:** If the parameter `--bootstrapping_rounds` is provided, it must be an integer between 0 and 3, with 0 having no effect.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/swgsrelate/usage) and the [parameter documentation](https://nf-co.re/swgsrelate/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/swgsrelate/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/swgsrelate/output).

## Credits

nf-core/swgsrelate was originally written by [Thomas Isensee](https://github.com/thomasisensee). This work was carried out as part of the [bwRSE4HPC](https://www.bwrse4hpc.de/) initiative, funded by the Baden-Württemberg Ministry of Science, Research and Arts, coordinated by the Scientific Software Center (SSC) at Heidelberg University and the Scientific Computing Center (SCC) at KIT.

We thank the following people for their extensive assistance in the development of this pipeline:

- [Gisela H. Kopp](https://github.com/GiselaHKopp)
- Till Dorendorf

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#swgsrelate` channel](https://nfcore.slack.com/channels/swgsrelate) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/swgsrelate for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
