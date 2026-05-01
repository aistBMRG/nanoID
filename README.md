
## nanoID: Reconstruction of Amplicon Sequence Variants and Species-Level Profiling Using Long Amplicon Reads


`nanoID` enables the reconstruction of Amplicon Sequence Variants (ASVs) from long‑read amplicon sequencing data and is applicable to medium‑ and high‑accuracy reads produced by Oxford Nanopore Technologies (ONT) and PacBio sequencing platforms. Further, `nanoID` enables robust species‑level profiling through greedy clustering of the reconstructed ASVs at 99% sequence identity, with subsequent quantitative inference of the species‑level operational taxonomic units (OTUs) using Emu.

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
wget [https://github.com/youruser/nanoid/releases/download/v0.1.0/nanoid-0.1.0.tar.gz](https://github.com/aistBMRG/nanoID/releases/download/v0.1.0/nanoid-0.1.0.tar.gz)
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
nanoid profile -i nanoid_condens -o nanoid_profile
```

**Defaults**

`nanoid condens`

```
General settings:
  --input_fastq INPUT_FASTQ, -i INPUT_FASTQ
                        Input FASTQ file containing raw Nanopore amplicon reads. (default: None)
  --output_directory OUTPUT_DIRECTORY, -o OUTPUT_DIRECTORY
                        Output directory for all generated files. (default: None)
  --randseed RANDSEED, -s RANDSEED
                        Random seed used for fastq splitting / subsampling. (default: 21336)
  --threads THREADS, -t THREADS
                        Number of threads. (default: 8)
  --keep_intermediates  Keep all intermediate files instead of deleting them. (default: False)

Cutadapt parameters:
  --cutadapt_f_primer CUTADAPT_F_PRIMER
                        Forward primer sequence. (default: AGRGTTYGATYHTGGCTCAG)
  --cutadapt_r_primer_rc CUTADAPT_R_PRIMER_RC
                        Reverse primer sequence (reverse-complemented). (default: AAGTCGTAACAAGGTARCCG)
  --cutadapt_f_primer_overlap CUTADAPT_F_PRIMER_OVERLAP
                        Min overlap for forward primer (None: use primer length). (default: None)
  --cutadapt_r_primer_overlap CUTADAPT_R_PRIMER_OVERLAP
                        Min overlap for reverse primer (None: use primer length). (default: None)
  --cutadapt_f_primer_errors CUTADAPT_F_PRIMER_ERRORS
                        Max mismatches for forward primer. (default: 2)
  --cutadapt_r_primer_errors CUTADAPT_R_PRIMER_ERRORS
                        Max mismatches for reverse primer. (default: 2)
  --cutadapt_mean_quality CUTADAPT_MEAN_QUALITY
                        Min mean read quality after trimming. (default: 20)
  --cutadapt_minimum_length CUTADAPT_MINIMUM_LENGTH
                        Min read length after primer trimming. (default: 1200)
  --cutadapt_maximum_length CUTADAPT_MAXIMUM_LENGTH
                        Max read length after primer trimming. (default: 1800)
  --cutadapt_maximum_expected_errors CUTADAPT_MAXIMUM_EXPECTED_ERRORS
                        Max expected errors after primer trimming. (default: 99)
  --skip_cutadapt       Skip primer trimming with Cutadapt. (default: False)

Splitting parameters:
  --disjoint_splits DISJOINT_SPLITS
                        Number of disjoint FASTQ splits. (default: 3)
  --max_reads_per_split MAX_READS_PER_SPLIT
                        Maximum number of reads per fastq split. (default: 40000)

Neighbor search parameters:
  --vsearch_id VSEARCH_ID
                        Min identity for near neighbor search. (default: 0.96)
  --vsearch_maxaccepts VSEARCH_MAXACCEPTS
                        Max number of near neighbors per read. (default: 10)

Condens parameters:
  --neighbors_n NEIGHBORS_N
                        Number of near neighbors per consensus. (default: 4)
  --kappa KAPPA         Dominance threshold for graph ascent. (default: 5.0)
  --ignore_singletons IGNORE_SINGLETONS
                        Ignore per-split singleton consensus sequences. (default: True)

Quant parameters:
  --min_identity MIN_IDENTITY
                        Min identity. (default: 0.98)
  --min_count MIN_COUNT
                        Min estimated read count for pruning. (default: 1.0)

Chimera check parameters:
  --abskew ABSKEW       abskew parameter for uchime3_denovo. (default: 4)
  --skip_uchime         Skip uchime3_denovo. (default: False)

```

`nanoid profile`

```
General settings:
  --threads THREADS, -t THREADS
                        Number of threads. (default: 8)
  --input_directory INPUT_DIRECTORY, -i INPUT_DIRECTORY
                        Input directory containing ASV files. (default: None)
  --output_directory OUTPUT_DIRECTORY, -o OUTPUT_DIRECTORY
                        Output directory for profiling results. (default: None)
  --fastq_pattern FASTQ_PATTERN
                        Filename pattern for FASTQ files. (default: _cutadapt_split.all.fastq)
  --fasta_pattern FASTA_PATTERN
                        Filename pattern for FASTA files with ASV sequences. (default: _uchime.fasta)

Import settings:
  --min_total_size_uniques MIN_TOTAL_SIZE_UNIQUES
                        Min total read count. (default: 0)
  --min_prev_uniques MIN_PREV_UNIQUES
                        Min prevalence. (default: 0)

Cluster settings:
  --cluster_id CLUSTER_ID
                        Identity threshold for clustering. (default: 0.99)
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

---

### Third-party attributions

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

### Citing nanoID

A preprint will be available shortly.

---

### Disclaimer

Portions of the codebase were developed with assistance from AI tools, specifically Claude AI (Sonnet 4.6) and Microsoft 365 Copilot (GPT‑5).

---

### License

Copyright (C) 2026 Dieter Tourlousse

This project is licensed under the GNU General Public License
version 3 or later (GPL-3.0-or-later).

See the LICENSE file for details.
