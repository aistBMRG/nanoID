
# nanoID: Reconstruction of amplicon sequence variants and species-level profiling using long amplicon reads

*nanoID* and its companying tool *nanoID_profile* reconstruct Amplicon Sequence Variants (AVSs) from long-read amplicon sequencing data. *nanoID* is applicable to medium-to-high accuracy reads generated using Oxford Nanopore Technologies (ONT) and PacBio.
---

## Installation

### Prerequisites

If `mamba` is not already installed:

```sh
conda install -n base -c conda-forge mamba
```

```sh
# Create environment and install dependencies
# Installation of Emu is optional
mamba create -n nanoid -c conda-forge -c bioconda \
    python=3.10 numpy scipy pandas biopython joblib rapidfuzz pyabpoa tqdm \
    vsearch cutadapt emu

# Activate environment
conda activate nanoid

# Download and install nanoID
wget https://github.com/youruser/nanoid/releases/download/v2.2.0/nanoid-2.2.0.tar.gz
pip install --no-build-isolation --no-deps nanoid-2.2.0.tar.gz 
rm nanoid-2.2.0.tar.gz
```

## Description of nanoID and nanoID_profile

## Quick start

```sh
nanoid -t 12 -i input.fastq -o nanoid
```
## Detailed usage

## Benchmarking data

## Publications

## License

Copyright (C) 2026 Dieter Tourlousse

This project is licensed under the GNU General Public License
version 3 or later (GPL-3.0-or-later).

See the LICENSE file for details.
