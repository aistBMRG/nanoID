
## nanoID: Reconstruction of Amplicon Sequence Variants and Species-Level Profiling Using Long Amplicon Reads

### What is nanoID?
`nanoID` is a new bioinformatics tool for accurate and sensitive reconstruction of amplicon sequence variants (ASVs) from long‑read amplicon sequencing data. The method is applicable to medium‑ and high‑accuracy reads generated using Oxford Nanopore Technologies (ONT) and PacBio sequencing platforms.
In addition to ASV reconstruction, `nanoID` enables robust species‑level profiling by greedy clustering of the reconstructed ASVs at 99% sequence identity, followed by quantification of the species‑level operational taxonomic units (OTUs) using Emu.

---

## 🚀 Installation

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
wget https://github.com/aistBMRG/nanoID/releases/download/v0.1.0/nanoid-0.1.0.tar.gz
pip install --no-build-isolation --no-deps nanoid-0.1.0.tar.gz
rm nanoid-0.1.0.tar.gz
```

##### Verify installation

After installation, confirm that `nanoid` is available:

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

## ⚡ Quick start

🖥️ **Basic usage**

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
nanoid profile -i nanoid_condens -o nanoid_profile
```

---

### 🔬 Schematic of `nanoid condens`

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

### 🔬 Schematic of `nanoid profile`

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

## Benchmarking/test data

In preparation.

---

## Third-party attributions

`nanoID` rely on several third‑party tools and we recommend citing the original publications of these tools.

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

## Citation

A manuscript describing nanoID is in preparation, and a preprint will be released shortly.  
In the meantime, please cite this repository when using nanoID.

---

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for a complete list of changes across releases.

---

## Disclaimer

Portions of the codebase were developed with assistance from AI tools, specifically Claude AI (Sonnet 4.6) and Microsoft 365 Copilot (GPT‑5).

---

## License

Copyright (C) 2026 Dieter Tourlousse

This project is licensed under the GNU General Public License
version 3 or later (GPL-3.0-or-later).

See the LICENSE file for details.
