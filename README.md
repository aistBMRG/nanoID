
## nanoID: Reconstruction of amplicon sequence variants and species-level profiling using long amplicon reads


`nanoID` reconstructs Amplicon Sequence Variants (AVSs) from long-read amplicon sequencing data and is applicable to medium-to-high accuracy reads generated using Oxford Nanopore Technologies (ONT) and PacBio. Its accompanying tool `nanoID_profile` combines and quantifies ASVs reconstructed from multiple samples, and also performs subsequent species-level clustering (99% identity) and quantification.

---

### Installation

The recommended way to install `nanoID` and its dependencies is via **conda** or **mamba**.

```sh
# Create a dedicated environment and install dependencies
# Note: Emu is optional and only required for species-level profiling
mamba create -n nanoid -c conda-forge -c bioconda \
    python=3.10 numpy scipy pandas biopython joblib rapidfuzz pyabpoa \
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

---

### Description of nanoID workflow

1. **Read preprocessing** (*optional*)  
   Raw reads are processed with **Cutadapt** for primer trimming and quality/lenghth filtering.

2. **Read partitioning**  
   Processed reads are divided into minimally overlapping partitions.  
   *By default, nanoID divides the reads into 2 non‑overlapping subsets.*

3. **Near‑neighbor identification**  
   For each partition, near‑neighbor reads are identified using **VSEARCH** (`--usearch_global`).

4. **Consensus sequence generation**  
   For each read within a partitioon, a consensus sequence is generated from its near‑neighbors using the
   abPOA partial‑order alignment algorithm.  
   *By default, nanoID uses 3 near-neighbors per read such that each consensus sequence is derived from 4 reads.*

5. **Denoising/condensing of consensus sequences**  
   Consensus sequences are denoised/condensed within each subset using nanoID’s *ConDens* algorithm, a greedy abundance‑based approach built on a phenomenological error model inspired by UNOISE3 (https://www.biorxiv.org/content/10.1101/081257v1).

6. **Recovery of ASVs**  
   Denoised/condensed consensus sequences consistently observed across partitions are retained as amplicon sequence variants (ASVs).

7. **Quantification of ASVs**  
   Abundances of each ASV are estimated using nanoID’s *Quant* algorithm, involving an expectation–maximization–like procedure to resolve read assignments.

8. **Chimera detection**  
   Putative chimeric ASVs are identified and removed using VSEARCH (`--uchime3_denovo`).

#### Output

The main output of `nanoID` is a fasta file of ASV sequences with USEARCH/VSEARCH-style size annotations.

```
>condenseq_1;size=3535
ACATAGATACAGTCTGATCGATCGTACGATCGATCGATCGATCGATCGATCGA
>condenseq_2;size=2148
ACATAGATACAGTCTGATCGATTGTACGATCGATCGATCGATCGATCGATCGA
>condenseq_3;size=987
ACATAGATACAGTCTGATCGATCGTACGATCGATCGATCGATCGATCGATCAA
>condenseq_4;size=412
ACATAGATACAGTCTGATCGATCGTACGATCGATCGATCGATCGATCGATTGA
```

### Description of nanoID_profile workflow

1. **Dereplication of ASVs accross samples** 
   Low abundance or infrequent ASVs are optionally removed.

2. **Clustering of ASVs in operational taxonomic units (OTUs)**  
   ASVs are greedily clustered in species-level OTUs VSEARCH's --cluster_smallmem at a global sequence identity of 99%, representing largely species-level OTUs.

3. **Quantification of OTU abundances** 
   The above generated OTUs are quantified using the processed reads with Emu (https://www.nature.com/articles/s41592-022-01520-4).

#### Output

The main output of `nanoID_profile` is a fasta file of OTU representatives and count table.

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
--cutadapt_r_primer_rc TGYACWCACCGCCCGTC
```

## Third-party attributions

`nanoID` and `nanoID_profile` rely on several third‑party tools and we recommend citing the original publications of these tools.

* Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016).
  VSEARCH: a versatile open source tool for metagenomics.
  PeerJ 4:e2584.
  doi:[10.7717/peerj.2584](https://doi.org/10.7717/peerj.2584)

* Curry KD, Wang Q, Nute MG, Tyshaieva A, Reeves E, Soriano S, Wu Q, Graeber E, Finzer P, Mendling W, Savidge T, Villapol S, Dilthey A, Treangen TJ (2022).
  Emu: species-level microbial community profiling of full-length 16S rRNA Oxford Nanopore sequencing data.
  Nat Methods 19(7):845-853.
  doi:[10.1038/s41592-022-01520-4](https://doi.org/10.1038/s41592-022-01520-4)

* Marcel M (2011).
  Cutadapt removes adapter sequences from high-throughput sequencing reads.
  EMBnet Journal 17(1):10-12.
  doi:[10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200)

## Citing nanoID

A preprint will be available shortly.

## License

Copyright (C) 2026 Dieter Tourlousse

This project is licensed under the GNU General Public License
version 3 or later (GPL-3.0-or-later).

See the LICENSE file for details.
