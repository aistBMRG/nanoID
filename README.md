
## nanoID: Reconstruction of amplicon sequence variants and species-level profiling using long amplicon reads


`nanoID` reconstructs Amplicon Sequence Variants (AVSs) from long-read amplicon sequencing data and is applicable to medium-to-high accuracy reads generated using Oxford Nanopore Technologies (ONT) and PacBio. Its accompanying tool `nanoID_profile` combines and quantifies ASVs reconstructed from multiple samples, and also performs subsequent species-level clustering (99% identity) and quantification.


### Installation

The recommended way to install `nanoID` and its dependencies is via **conda** or **mamba**.

```sh
# Create a dedicated environment and install dependencies
# Note: Emu is optional and only required for species-level profiling
mamba create -n nanoid -c conda-forge -c bioconda \
    python=3.10 numpy scipy pandas biopython joblib rapidfuzz pyabpoa tqdm \
    vsearch cutadapt emu

# Activate the environment
mamba activate nanoid

# Download and install nanoID
wget https://github.com/youruser/nanoid/releases/download/v2.2.0/nanoid-2.2.0.tar.gz
pip install --no-build-isolation --no-deps nanoid-2.2.0.tar.gz
rm nanoid-2.2.0.tar.gz
```

##### Verify installation

After installation, confirm that `nanoID` and `nanoID_profile` are available by checking their help messages:

```sh
nanoid -h
nanoid_profile -h
```

### Description of nanoID

`nanoID` consists of the following steps:
1. *(Optional)* Processing of reads using `Cutadapt` for primer trimming and filtering.
2. Partitioning of processed reads into minimally overlapping subsets. By default, `nanoID` partitions the reads into 3 non-overlapping subsets. 
3. For each subset, identifying near‑neighbors using `VSEARCH` (--usearch_global command).
4. For each subset, generating consensus sequences from each read using the abPOA partial order alignment algorithm. By default, `nanoID` uses 3 near-neigbors per read. 
5. For each subset, condensing (denoising) of the consensus sequences using `nanoid`'s ConDens algorithm.
6. Recovering of condensed consensus sequences consistently observed across partitions as *trusted* ASVs.
7. Quantifying the ASVs using using `nanoid`'s Quant algorithm.
8. Detecting chimeras using `VSEARCH` (--uchime3_denovo command).


## Quick start

```sh
nanoid -t 12 -i input.fastq -o nanoid
```

## Detailed usage

## Benchmarking data

nanoID was tested using publically available sequencing data for Zymo's mock community, with default settings except for Cutadapt.

```sh
nanoid -t 12 \
--cutadapt_f_primer    AGRGTTYGATYMTGGCTCAG \
--cutadapt_r_primer_rc TGYACWCACCGCCCGTC \
--cutadapt_minimum_length 750 \
--cutadapt_maximum_length 1750
```

## Publications

## License

Copyright (C) 2026 Dieter Tourlousse

This project is licensed under the GNU General Public License
version 3 or later (GPL-3.0-or-later).

See the LICENSE file for details.
