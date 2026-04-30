
## nanoID: Reconstruction of Amplicon Sequence Variants and Species-Level Profiling Using Long Amplicon Reads


`nanoID` enables the reconstruction of Amplicon Sequence Variants (ASVs) from long‑read amplicon sequencing data and is applicable to medium‑ and high‑accuracy reads produced by Oxford Nanopore Technologies (ONT) and PacBio sequencing platforms. Further, `nanoID` enables robust species‑level profiling through greedy clustering of ASVs at 99% sequence identity, with subsequent quantitative inference of species‑level operational taxonomic units (OTUs) using Emu.

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
wget https://github.com/youruser/nanoid/releases/download/v0.1.0/nanoid-0.1.0.tar.gz
pip install --no-build-isolation --no-deps nanoid-0.1.0.tar.gz
rm nanoid-0.1.0.tar.gz
```

##### Verify installation

After installation, confirm that `nanoID` and `nanoID_profile` are available by checking their help messages:

```sh
nanoid -h
```

```
usage: nanoid [-h] [-v] <command> ...

nanoID: Long‑read amplicon denoising and species profiling

options:
  -h, --help     show this help message and exit
  -v, --version  Show program version and exit.

Subcommands:
  <command>
    condens      Workflow for recovering ASVs.
    profile      Workflow for species-level profiling based on recovered ASVs.
```

---

### Quick start

**Basic usage**

Run `nanoid condens` on a single sample:
```sh
nanoid condens -i input.fastq -o nanoid_condens
```

For downstream joint profiling across multiple samples, write all `nanoid` outputs to a single directory.

```sh
mkdir nanoid

## Using a shell loop
for fq in *.fastq; do
    nanoid condens -i "$fq" -o nanoid_condens
done

## Using GNU Parallel
parallel -j 4 \
    "nanoid condens -i {} -o nanoid_condens" \
    ::: *.fastq
```

Once all samples have been processed, run `nanoid profile` on the combined output directory.
```sh
nanoid profile -t 12 -i nanoid_condens -o nanoid_profile
```

---

### Schematic of `nanoid condens`

```
reads
↓
[cutadapt, optional] ─ primer trimming, length/quality filtering, read orientation (--revcomp flag)
↓
[N disjoint_splits]
  ├─ split 1 ─▶ [vsearch, all vs. all] ─▶ [abPOA consensus, neighbors_n reads] ─▶ conseqs₁
  ├─ split 2 ─▶ [vsearch, all vs. all] ─▶ [abPOA consensus, neighbors_n reads] ─▶ conseqs₂
  └─ split N ─▶ [vsearch, all vs. all] ─▶ [abPOA consensus, neighbors_n reads] ─▶ conseqsₙ
↓
┌─────────────────────────────────────────────┐
  CROSS‑SPLIT CONSENSUS GRAPH
  nodes = intersection(conseqs₁ … conseqsₙ)
  edges = union(shared‑neighbor edges)
└─────────────────────────────────────────────┘
↓
[graph ascent, κ threshold]
  cross-split constrained graph ascent
↓
[graph peaks] → condensed conseqs
↓
[EM quantification] → read assignment
↓
[uchime3_denovo] → chimera removal
↓
Amplicon Sequence Variants (ASVs)
```
##### Output

* Fasta file `{basename}_nanoid_uchime.fasta` of ASV sequences with USEARCH/VSEARCH-style size annotations.
* Log file `{basename}_nanoid_condens.log`

### Schematic of `nanoid profile`

```
ASVs
↓
[dereplication across samples]
  └─ optional: total count / prevalence filtering
↓
[VSEARCH clustering]
  └─ greedy clustering into OTUs [nanOTUs]
     (--cluster_smallmem, 99% identity)
↓
[nanOTUs]
↓
[Emu quantification]
↓
nanOTU abundance profiles
```
##### Output

* Fasta file `{basename}_nanoid_nanotus.fasta` of OTU sequences.
* Count table `{basename}_nanoid_nanotus_counts.tsv` with abundances estimated using Emu against the OTU sequences.
* Log file `{basename}_nanoid_profile.log`

---

### Benchmarking/test data

We applied `nanoID` to ONT sequencing data (R10.4.1 flowcell, basecalling with [Dorado](https://github.com/nanoporetech/dorado) in super accuracy mode, model sup@5.2.0) for the ZymoBIOMICS Microbial Community DNA Standard (catalog no. D6305), as described by [Riisgaard-Jensen et al.](https://www.biorxiv.org/content/10.64898/2026.02.26.708165v1). Data for replicate 1 was obtained from SRA under accession number [SRR36567714](https://www.ncbi.nlm.nih.gov/sra/?term=SRR36567714). 

```sh
prefetch SRR36567714
fasterq-dump SRR36567714
```

```sh
nanoid condens -t 8 \
-i SRR36567714.fastq \
-o SRR36567714_nanoid \
--max_reads_per_split 20000 \
--cutadapt_f_primer    AGRGTTYGATYMTGGCTCAG \
--cutadapt_r_primer_rc TGYACWCACCGCCCGTC
```

This yielded <u><strong>27 non-chimeric ASVs</strong></u> 
that perfectly matched:
* [ZymoBIOMICS.STD.refseq.v2.fasta](https://s3.amazonaws.com/zymo-files/BioPool/ZymoBIOMICS.STD.refseq.v2.zip) reference sequences (N=17).
* Sequences from same species in NCBI (N=10), suggesting that they represented genuine sequences.

```sh
nanoid_profile -t 8 \
-i SRR36567714_nanoid \
-o SRR36567714_nanoid_profile \
--cluster_id 0.99
```
This yielded **9** OTUs.

---

### Third-party attributions

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

---

### Citing nanoID

A preprint will be available shortly.

---

### License

Copyright (C) 2026 Dieter Tourlousse

This project is licensed under the GNU General Public License
version 3 or later (GPL-3.0-or-later).

See the LICENSE file for details.
