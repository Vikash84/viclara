# ViClaRA
![](http://img.shields.io/badge/nextflow-v20.04.1-lightgreen)
![](https://img.shields.io/badge/uses-docker-blue.svg)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/mlconjug/badges/installer/conda.svg)](https://conda.anaconda.org/conda-forge)
![](https://img.shields.io/badge/licence-MIT-lightgrey.svg)

[![Twitter Follow](https://img.shields.io/twitter/follow/john_juma.svg?style=social)](https://twitter.com/john_juma)

We present [ViClaRA], a nextflow pipeline to classify reads and perform reference guided assembly on viral metagenomics reads from Illumina platform. Metagenomics and high throughput genomics methods allow us to identify potential reservoir of zoonotic pathogens and discover new variants. Metagenomics is a powerful approach for the broad identification of pathogens in clinical samples.

## Installation

You only need [Nextflow](https://nf-co.re/usage/installation) (version 20.+) [Docker](https://docs.docker.com/engine/installation/) and [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) installed to run the pipeline. All dependencies will be pulled automatically. 

Either run ViClaRA by cloning this repository:
```bash
git clone https://github.com/ajodeh-juma/viclara.git
cd viclara
nextflow run viclara.nf --help
```

or let Nextflow do the pull
```bash
nextflow pull ajodeh-juma/viclara
```