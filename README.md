
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

nanoID – long‑read amplicon denoising and profiling

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

Run nanoid on a single sample:
```sh
nanoid -t 12 -i input.fastq -o nanoid
```
This generates per‑sample results in the specified output directory.

For downstream joint profiling across multiple samples, write all `nanoid` outputs to a single directory.

```sh
mkdir nanoid

## Using a shell loop
for fq in *.fastq; do
    nanoid -t 12 -i "$fq" -o nanoid
done

## Using GNU Parallel
parallel -j 4 \
    "nanoid -t 12 -i {} -o nanoid" \
    ::: *.fastq
```

Once all samples have been processed, run `nanoid_profile` on the combined output directory.
```sh
nanoid_profile -t 12 -i nanoid -o nanoid_profile
```

---

#### Description of `nanoid condens` workflow

1. **Read preprocessing** (*optional*)  
   Raw reads are processed with [Cutadapt](https://github.com/marcelm/cutadapt) for primer trimming and quality/lenghth filtering.

2. **Read partitioning**  
   Reads are divided into disjoint partitions.  
   > By default, nanoID divides the reads into 3 partitions.

3. **Near‑neighbor identification**  
   For each partition, near‑neighbor reads are identified using [VSEARCH](https://github.com/torognes/vsearch) (`--usearch_global`).

4. **Generation of consensus sequences (conseqs)**  
   For each read within a partition, a consensus sequence is generated from its near‑neighbors using the
   [abPOA](https://github.com/yangao07/abpoa) partial‑order alignment algorithm.  
   > By default, nanoID uses 4 near-neighbors per read such that each consensus is derived from 5 reads.

5. **Denoising/condensing of conseqs**  
  Conseqs are denoised/condensed using nanoID’s *Condens* algorithm. Briefly, *Condens* dereplicates consensus sequences and constructs a shared‑neighbor graph for each split, where nodes represent consensus sequences and edges reflect shared read support. A cross‑split consensus graph is then derived, retaining (i) only nodes observed in all splits and (ii) all edges observed in at least one split. This consensus graph is subsequently subjected to κ‑ascent to identify a set of well‑supported peaks (condensed conseqs), which are retained as amplicon sequence variants (ASVs). To prevent the loss of closely related ASVs that share neighbors, ascent is constrained by a minimum abundance ratio between connected nodes, defined by the parameter κ.

6. **Quantification of ASVs**  
   Abundances of the ASVs are estimated using nanoID’s *Quant* algorithm, an expectation–maximization–like procedure to resolve ambiguous read assignments.

8. **Chimera detection**  
   Putative chimeric ASVs are identified and removed using VSEARCH (`--uchime3_denovo`).

##### Output

* Fasta file `{basename}_nanoid_nonchimeras.fasta` of ASV sequences with USEARCH/VSEARCH-style size annotations.
* Log file `{basename}_nanoid.log`

#### Description of `nanoid profile` workflow

1. **Dereplication of ASVs accross samples**  
   Low abundance or infrequent ASVs are optionally removed.

3. **Clustering of ASVs in operational taxonomic units (OTUs)**  
   ASVs are greedily clustered in OTUs using VSEARCH's `--cluster_smallmem` at a global sequence identity of 99%, representing largely species-level clusters.

4. **Quantification of OTU abundances**  
   OTUs are quantified using the processed reads with [Emu](https://www.nature.com/articles/s41592-022-01520-4).

##### Output

* Fasta file `{basename}_nanoid_centroids.fasta` of ASV sequences with USEARCH/VSEARCH-style size annotations.
* Count table `{basename}nanoid_centroids_profile_emu.tsv`
* Log file `{basename}_nanoid_profile.log`

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
--disjoint_partitions 3 \
--max_reads_per_disjoint_partition 50000 \
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
