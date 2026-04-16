
## nanoID: Reconstruction of amplicon sequence variants and species-level profiling using long amplicon reads


`nanoID` reconstructs Amplicon Sequence Variants (AVSs) from long-read amplicon sequencing data and is applicable to medium-to-high accuracy reads generated using Oxford Nanopore Technologies (ONT) and PacBio. Its accompanying tool `nanoID_profile` enables species-level profiling by greedy clustering of the ASVs at 99% sequence identity, and robust quantification of the species-level operational taxonomic units (OTUs) using Emu (Expectation–Maximization-based taxonomic Unit profiling).

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
nanoid_profile -h
```

---

### Quick start

**Basic usage**
```sh
nanoid -t 12 -i input.fastq -o nanoid
```

For combined analysis of multiple samples, the output of `nanoid` can be written to a single output directory.

```sh
mkdir nanoid

for fq in *.fastq; do
    nanoid -t 12 -i "$fq" -o nanoid
done

## with GNU parallel

parallel -j 4 \
    "nanoid -t 12 -i {} -o nanoid" \
    ::: *.fastq
```

```sh
nanoid_profile -t 12 -i nanoid -o nanoid_profile
```

---

#### Description of nanoID workflow

1. **Read preprocessing** (*optional*)  
   Raw reads are processed with **Cutadapt** for primer trimming and quality/lenghth filtering.

2. **Read partitioning**  
   Reads are divided into minimally overlapping partitions.  
   *By default, nanoID divides the reads into 2 non‑overlapping subsets.*

3. **Near‑neighbor identification**  
   For each partition, near‑neighbor reads are identified using VSEARCH (`--usearch_global`).

4. **Consensus sequence generation**  
   For each read within a partition, a consensus sequence is generated from its near‑neighbors using the
   abPOA partial‑order alignment algorithm.  
   *By default, nanoID uses 3 near-neighbors per read such that each consensus is derived from 4 reads.*

5. **Denoising/condensing of consensus sequences**  
   Consensus sequences are denoised/condensed within each partition using nanoID’s *ConDens* algorithm, a greedy abundance‑based approach built on a phenomenological error model inspired by [UNOISE3](https://www.biorxiv.org/content/10.1101/081257v1).

6. **Recovery of ASVs**  
   Denoised/condensed consensus sequences consistently observed across partitions are retained as amplicon sequence variants (ASVs).

7. **Quantification of ASVs**  
   ASV abundances are estimated using nanoID’s *Quant* algorithm,an expectation–maximization–like procedure to resolve read assignments.

8. **Chimera detection**  
   Putative chimeric ASVs are identified and removed using VSEARCH (`--uchime3_denovo`).

##### Output

* Fasta file `{basename}_nanoid_nonchimeras.fasta` of ASV sequences with USEARCH/VSEARCH-style size annotations.
* Log file `{basename}_nanoid.log`

#### Description of nanoID_profile workflow

1. **Dereplication of ASVs accross samples**  
   Low abundance or infrequent ASVs are optionally removed.

3. **Clustering of ASVs in operational taxonomic units (OTUs)**  
   ASVs are greedily clustered in OTUs using VSEARCH's --cluster_smallmem at a global sequence identity of 99%, representing largely species-level clusters.

4. **Quantification of OTU abundances**  
   OTUs are quantified using the processed reads with [Emu](https://www.nature.com/articles/s41592-022-01520-4).

##### Output

* Fasta file `{basename}_nanoid_centroids.fasta` of ASV sequences with USEARCH/VSEARCH-style size annotations.
* Count table `{basename}nanoid_centroids_profile_emu.tsv`
* Log file `{basename}_nanoid_profile.log`

---

## ConDens: condensation algorithm

`ConDens` operates on dereplicated, abundance‑sorted consensus sequences (*children*) {c<sub>1</sub>, c<sub>2</sub>, …, c<sub>n</sub>} with abundances {a<sub>1</sub>, a<sub>2</sub>, …, a<sub>n</sub>}, where a<sub>1</sub> ≥ a<sub>2</sub> ≥ … ≥ a<sub>n</sub>.

The goal is to generate a reduced set of denoised / condensed consensus
sequences (*parents*) {p<sub>1</sub>, p<sub>2</sub>, …, p<sub>k</sub>}, with abundances {A<sub>1</sub>, A<sub>2</sub>, …, A<sub>n</sub>}, by greedily assigning children to parents based on abundance and sequence divergence. The algorithm is inspired by UNOISE3.

Sequence divergence between two sequences is measured as the lenght-normalized Levenshtein edit distance:
> d(s<sub>1</sub>, s<sub>2</sub>) = D(s<sub>1</sub>, s<sub>2</sub>) / max(|s<sub>1</sub>|, |s<sub>2</sub>|)


##### Initialization

The most abundant child becomes the first parent:

> s<sub>1</sub> → p<sub>1</sub>

##### Parent–child scoring function

For each candidate parent p<sub>i</sub> of a given child s<sub>j</sub>,
ConDens computes score 𝒮(p<sub>i</sub>, s<sub>j</sub>) as:

> 𝒮(p<sub>i</sub>, s<sub>j</sub>) = A<sub>i</sub> · exp(−α · d(p<sub>i</sub>, s<sub>j</sub>)<sup>γ</sup> − log β)

where:
- A<sub>i</sub> is the abundance of parent p<sub>i</sub>,
- β enforces a baseline parent‑to‑child abundance ratio (default: 8),
- α controls the strength of the sequence divergence penalty (default: 0.5),
- γ controls the non‑linearity of the divergence penalty (default: 1.5).

##### Greedy assignment rule

For a given child s<sub>j</sub>, the best parent is selected as:

> p*<sub>j</sub> = arg max<sub>p<sub>i</sub></sub> 𝒮(p<sub>i</sub>, s<sub>j</sub>)

The child is assigned to this parent if:

> 𝒮(p*<sub>j</sub>, s<sub>j</sub>) ≥ a<sub>j</sub>

where a<sub>j</sub> is the abundance of child s<sub>j</sub>.

If no parent satisfies this condition, the child becomes a new parent:

> s<sub>j</sub> → p<sub>k+1</sub>

##### Parent reinforcement (optional)

If parent reinforcement is enabled, accepted assignments update the parent abundance:

> A<sub>p*</sub> ← A<sub>p*</sub> + a<sub>j</sub>

---

### Benchmarking/test data

We applied `nanoID` to ONT sequencing data (R10.4.1 flowcell, basecalling with [Dorado](https://github.com/nanoporetech/dorado) in super accuracy mode, model sup@5.2.0) for the ZymoBIOMICS Microbial Community DNA Standard (catalog no. D6305), as described by [Riisgaard-Jensen et al.](https://www.biorxiv.org/content/10.64898/2026.02.26.708165v1). Data for replicate 1 was obtained from SRA under accession number [SRR36567714](https://www.ncbi.nlm.nih.gov/sra/?term=SRR36567714). 

```sh
prefetch SRR36567714
fasterq-dump SRR36567714
```

```sh
# nanoID with 2x 50,000 randomly selected reads, after processing using Cutadapt #
nanoid -t 8 \
-i SRR36567714.fastq \
-o SRR36567714_nanoid \
--number_fastq_splits 2 \
--reads_per_split 50000 \
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
This yielded **9** OTUs, including two ASVs affilitated for .

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
