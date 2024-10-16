# LZ-ANI

[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/lz-ani/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/LZ-ANI/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/lz-ani.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/lz-ani)
[![GitHub Actions CI](../../workflows/GitHub%20Actions%20CI/badge.svg)](../../actions/workflows/main.yml)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

LZ-ANI is a fast and memory-efficient tool for determining average nucleotide identity (ANI) in large sets of genomic sequences. The tool uses Lempel-Ziv parsing, achieving high sensitivity in discerning matched and mismatched nucleotides, thereby enabling accurate determination of ANI. Its efficiency stems from a simplified indel handling model, making LZ-ANI magnitudes faster than alignment-based tools (e.g., BLASTn, MegaBLAST) while maintaining accuracy comparable to the most sensitive BLASTn searches. 

LZ-ANI is a key component of [Vclust](https://github.com/refresh-bio/vclust), a tool for comparing and clustering virus genomes. Although optimized for virus genomes, LZ-ANI's parameters can be customized for longer genomes, such as those of bacteria and archaea.

#### Sequence similarity measures

LZ-ANI offers six similarity measures between two genomic sequences:

- **ANI**: The number of identical bases across local alignments divided by the total length of the alignments.
- **Global ANI (gANI)**: The number of identical bases across local alignments divided by the length of the query/reference genome.
- **Total ANI (tANI)**: The number of identical bases between query-reference and referece-query genomes divided by the sum length of both genomes.
- **Coverage (alignment fraction)**: The proportion of the query/reference sequence aligned with the reference/query sequence.
- **Number of local alignments**: The count of individual alignments found between the sequences.
- **Ratio between query and reference genome lengths**: A measure comparing the lengths of the two genomes.


## Installation

### Prebuilt binaries
The easiest way to get started is by downloading the prebuilt binaries from the Releases page. Select your platform and download the tool.

### Building from source
Alternatively, you can clone the repository and build LZ-ANI from source. You will need a modern C++17 compiler (e.g., g++10 or newer) and NASM, which can often be installed via a system package manager:

```bash
sudo apt install nasm                                     # Ubuntu
```

Then, clone the repository and build the LZ-ANI binary:

```bash
git clone --recurse-submodules https://github.com/refresh-bio/lz-ani
cd lz-ani && make -j
```

## Quick start

Perform pairwise alignment among all genome pairs within a single multifasta file:

```bash
./lz-ani all2all --in-fasta example/multifasta.fna --out ani.tsv
```

Perform pairwise alignment for all genome files in a directory:

```bash
./lz-ani all2all --in-dir example/fna --out ani.tsv
```

## Usage

```bash
lz-ani <mode> [options]
```

#### Modes
Currently, LZ-ANI operates exclusively in the all2all mode, where sequence similarities are computed for all pairs of input sequences. It supports filtering the list of pairs, with future updates to introduce additional modes.

#### Input files specification:
* `--in-fasta <file_name>` &mdash; FASTA file (for multisample-fasta mode)
* `--in-txt <file_name>` &mdash; text file with FASTA file names
* `--in-dir <path>` &mdash; directory with FASTA files
* `--multisample-fasta <bool>` &mdash; multi sample FASTA input (default: true)
* `--flt-kmerdb <fn> <float>` &mdash; filtering file (kmer-db output) and threshold

#### Output files specification:
* `-o, --out <file_name>` &mdash; output file name
* `--out-ids <file_name>` &mdash; output file name for ids file (optional)
* `--out-in-percent <bool>` &mdash; output in percent (default: false)
* `--out-type <type>` &mdash; one of:
  * `tsv` &mdash; two tsv files with: results defined by `--out-format` and sequence ids (default)
  * `single-txt` &mdash; combined results in single txt file
* `--out-format <type>` &mdash; comma-separated list of values:
  *  `id1`, `id2`, `idx1`, `idx2`, `len1`, `len2`, `tani`, `gani`, `ani`, `cov`, `len_ratio`, `nt_match`, `nt_mismatch`, `num_alns`
  * you can include also meta-names:
    * `complete=idx1,idx2,id1,id2,tani,gani,ani,cov,num_alns,len_ratio,len1,len2,nt_match,nt_mismatch`
    * `lite=idx1,idx2,tani,gani,ani,cov,num_alns,len_ratio`
    * `standard=idx1,idx2,id1,id2,tani,gani,ani,cov,num_alns,len_ratio`
    * `(default: standard)`
* `--out-alignment <file_name>` &mdash; output file name for alignments (optional)
* `--out-filter <par> <float>` &mdash; store only results with `<par>` (can be: `tani`, `gani`, `ani`, `cov`) at least `<float>`; can be used multiple times

#### LZ-parsing options:
* `-a, --mal <int>` &mdash; min. anchor length (default: 11)
* `-s, --msl <int>` &mdash; min. seed length (default: 7)
* `-r, --mrd <int>` &mdash; max. dist. between approx. matches in reference (default: 40)
* `-q, --mqd <int>` &mdash; max. dist. between approx. matches in query (default: 40)
* `-g, --reg <int>` &mdash; min. considered region length (default: 35)
* `--aw <int>` &mdash; approx. window length (default: 15)
* `--am <int>` &mdash; max. no. of mismatches in approx. window (default: 7)
* `--ar <int>` &mdash; min. length of run ending approx. extension (default: 3)

#### Other options:
* `-t, --threads <int>` &mdash; no of threads; 0 means auto-detect (default: 0)
* `-V, --verbose <int>` &mdash; verbosity level (default: 1)
  


## Input
LZ-ANI accepts a single FASTA file containing genomic sequences ([example](./example/multifasta.fna)) or a directory of FASTA files (one genome per file) ([example](./example/fna)). The input file(s) can be gzipped.

```bash
# Perform pairwise alignment among all genome pairs within a single multifasta file.
./lz-ani all2all --in-fasta example/multifasta.fna --out ani.tsv

# Perform pairwise alignment for all genome files in a directory.
./lz-ani all2all --in-dir example/fna --out ani.tsv
```



## Output

LZ-ANI creates two TSV files: one contains ANI values for genome pairs, and the other lists genome identifiers sorted by decreasing sequence length. For example, the following command will create two TSV files: [ani.tsv](./example/output/ani.tsv) and [ani.ids.tsv](./example/output.ani.ids.tsv):


```bash
./lz-ani all2all --in-fasta example/multifasta.fna --out ani.tsv
```

For brevity, only the first 12 lines of output are shown:

```
qidx  ridx  query reference   tani  gani  ani   qcov  rcov  num_alns len_ratio
9  8  NC_010807.alt3 NC_010807.alt2 0.972839 0.960192 0.986657 0.973177 0.997608 60 0.9836
8  9  NC_010807.alt2 NC_010807.alt3 0.972839 0.985279 0.987642 0.997608 0.973177 67 0.9836
10 8  NC_010807.alt1 NC_010807.alt2 0.987250 0.987041 0.987117 0.999923 0.999901 34 0.9571
8  10 NC_010807.alt2 NC_010807.alt1 0.987250 0.987449 0.987547 0.999901 0.999923 36 0.9571
11 8  NC_010807.ref  NC_010807.alt2 0.989807 0.989540 0.989617 0.999923 1.000000 14 0.9571
8  11 NC_010807.alt2 NC_010807.ref  0.989807 0.990063 0.990063 1.000000 0.999923 14 0.9571
10 9  NC_010807.alt1 NC_010807.alt3 0.979963 0.993250 0.994557 0.998686 0.972575 71 0.9730
9  10 NC_010807.alt3 NC_010807.alt1 0.979963 0.967035 0.994304 0.972575 0.998686 70 0.9730
11 9  NC_010807.ref  NC_010807.alt3 0.983839 0.997166 0.997217 0.999948 0.974230 52 0.9730
9  11 NC_010807.alt3 NC_010807.ref  0.983839 0.970871 0.996552 0.974230 0.999948 52 0.9730
11 10 NC_010807.ref  NC_010807.alt1 0.997462 0.997475 0.997475 1.000000 1.000000 23 1.0000
10 11 NC_010807.alt1 NC_010807.ref  0.997462 0.997449 0.997449 1.000000 1.000000 23 1.0000
```

### Output format

The `--out-format` provides three output views: `standard`, `lite`, and `complete`.

| Column | Standard | Lite | Complete | Description |
| --- | :---: |:---: | :---: | --- |
| qidx | + | + | +  | Index of query sequence |
| ridx | + | + | +  | Index of reference sequence |
| query | + | - | +  | Identifier (name) of query sequence |
| reference | + | - | +  | Identifier (name) of reference sequence |
| tani | + | + | +  | total ANI [0-1] |
| gani | + | + | +  | global ANI [0-1] |
| ani | + | + | +  | ANI [0-1] |
| qcov | + | + | +  | Query coverage (aligned fraction) [0-1] |
| rcov | + | + | +  | Reference coverage (aligned fraction) [0-1] |
| num_alns | + | + | +  | Number of local alignments |
| len_ratio | + | + | +  | Length ratio between shorter and longer sequence [0-1] |
| qlen | - | - | +  | Query sequence length |
| rlen | - | - | +  | Reference sequence length |
| nt_match | - | - | +  | Number of matching nucleotides across alignments |
| nt_mismatch | - | - | +  | Number of mismatching nucleotides across alignments |


In addition, the `--out-format` option permits formatting arbitrary fields from the LZ-ANI tab-separated-value (TSV) format: 

```bash
./lz-ani all2all --in-fasta example/multifasta.fna --out ani.tsv --out-format query,reference,ani,qcov,rcov
```

```
query reference ani qcov  rcov
NC_010807.alt2  NC_025457.alt2  0.572519  0.0646036 0.0387601
NC_025457.alt2  NC_010807.alt2  0.570567  0.0387601 0.0646036
NC_010807.alt3  NC_025457.alt2  0.586745  0.0514402 0.0354560
NC_025457.alt2  NC_010807.alt3  0.565714  0.0354560 0.0514402
NC_010807.alt1  NC_025457.alt2  0.577825  0.0604148 0.0394770
NC_025457.alt2  NC_010807.alt1  0.568496  0.0394770 0.0604148
NC_010807.ref NC_025457.alt2  0.57375 0.0618318 0.0395705
NC_025457.alt2  NC_010807.ref 0.567546  0.0395705 0.0618318
NC_005091.alt1  NC_005091.alt2  0.937913  0.996571  0.996907
NC_005091.alt2  NC_005091.alt1  0.940487  0.996907  0.996571
NC_005091.ref NC_005091.alt2  0.964911  0.999495  0.999859
NC_005091.alt2  NC_005091.ref 0.968125  0.999859  0.999495
NC_002486.alt NC_005091.alt2  0.558574  0.0129065 0.00871326
...
```


### Output filtering

The `--out-filter` option allows you to filter the output by setting minimum similarity thresholds, enabling you to report only those genome pairs that meet the specified criteria, thus significantly reducing the output TSV file size. For example, the following command outputs only genome pairs with ANI ≥ 0.95 and query coverage ≥ 0.85:

```bash
./lz-ani all2all --in-fasta example/multifasta.fna --out ani.tsv --out-filter ani 0.95 --out-filter qcov 0.85
```

```
qidx  ridx  query reference tani  gani  ani qcov  num_alns  len_ratio
7 6 NC_025457.alt1  NC_025457.ref 0.809496  0.845785  0.985613  0.858131  123 0.9628
9 8 NC_010807.alt3  NC_010807.alt2  0.972839  0.960192  0.986657  0.973177  60  0.9836
8 9 NC_010807.alt2  NC_010807.alt3  0.972839  0.985279  0.987642  0.997608  67  0.9836
10  8 NC_010807.alt1  NC_010807.alt2  0.987250  0.987041  0.987117  0.999923  34  0.9571
8 10  NC_010807.alt2  NC_010807.alt1  0.987250  0.987449  0.987547  0.999901  36  0.9571
11  8 NC_010807.ref NC_010807.alt2  0.989807  0.989540  0.989617  0.999923  14  0.9571
8 11  NC_010807.alt2  NC_010807.ref 0.989807  0.990063  0.990063  1 14  0.9571
10  9 NC_010807.alt1  NC_010807.alt3  0.979963  0.993250  0.994557  0.998686  71  0.9730
9 10  NC_010807.alt3  NC_010807.alt1  0.979963  0.967035  0.994304  0.972575  70  0.9730
11  9 NC_010807.ref NC_010807.alt3  0.983839  0.997166  0.997217  0.999948  52  0.9730
9 11  NC_010807.alt3  NC_010807.ref 0.983839  0.970871  0.996552  0.974230  52  0.9730
11  10  NC_010807.ref NC_010807.alt1  0.997462  0.997475  0.997475  1 23  1
10  11  NC_010807.alt1  NC_010807.ref 0.997462  0.997449  0.997449  1 23  1
...
```

## Alignments between selected genome pairs

By default, LZ-ANI performs all-to-all pairwise sequence alignments. Alternatively, it can align only specific sequence pairs provided in a user-defined text filter file. This file ([example](./example/fltr.txt)), which can be created using [Kmer-db](https://github.com/refresh-bio/kmer-db) or manually, contains the pairs for alignment. To align only those genome pairs with similarity values of 0.9 or greater, use the following command:

```bash
./lz-ani all2all --in-fasta example/multifasta.fna --out ani.tsv --flt-kmerdb example/fltr.txt 0.9
```

### Create a filter file using Kmer-db

1. Build a database from all *k*-mers contained in a FASTA file of genome sequences.

```bash
# Create a text file containing a path to input FASTA file.
ls example/multifasta.fna > ids.txt
# Build a database from 25-mers.
kmer-db build -multisample-fasta -k 25 ids.txt multifasta.kdb
```

2. Identify number of common *k*-mers between all sequence pairs in the database.

```bash
# Output sequence pairs with more than 10 common 25-mers, outputting in a sparse matrix.
kmer-db all2all -sparse -above 10 multifasta.kdb all2all.txt
```

3. Compute sequence identity (relative to the shorter sequence) for genome pairs with the output from the previous step, generating a filter file for LZ-ANI:

```bash
# Calculate sequence identity and create the filter file.
kmer-db distance ani-shorter -sparse -above 0.7 all2all.txt
mv all2all.txt fltr.txt
```

### Alignments

LZ-ANI can output alignment details in a separate TSV file. This output format is similar to the BLASTn tabular output and includes information on each local alignment between two genomes, such as the coordinates in both the query and reference sequences, strand orientation, the number of matched and mismatched nucleotides, and the percentage of sequence identity.

```bash
./lz-ani all2all --in-fasta example/multifasta.fna --out ani.tsv --out-alignment ani.aln.tsv
```

Sample output:

```
query reference   pident   alnlen   qstart   qend  rstart   rend  nt_match nt_mismatch
NC_025457.alt2 NC_025457.ref  89.2893  999   22119 23117 14207 15163 892   107
NC_025457.alt2 NC_025457.ref  89.8305  826   3373  4198  2202  3020  742   84
NC_025457.alt2 NC_025457.ref  91.0804  796   41697 42492 27680 28475 725   71
NC_025457.alt2 NC_025457.ref  87.2483  745   38039 38783 24969 25688 650   95
NC_025457.alt2 NC_025457.ref  89.8860  702   7269  7970  5077  5778  631   71
NC_025457.alt2 NC_025457.ref  93.2081  692   62572 63263 41329 42020 645   47
NC_025457.alt2 NC_025457.ref  90.9565  575   31121 31695 20438 21003 523   52
NC_025457.alt2 NC_025457.ref  90.6195  565   11476 12040 7999  8563  512   53
NC_025457.alt2 NC_025457.ref  91.6211  549   10905 11453 7455  8003  503   46
NC_025457.alt2 NC_025457.ref  86.7041  534   29624 30157 19067 19586 463   71
NC_025457.alt2 NC_025457.ref  93.5673  513   10149 10661 6915  7427  480   33
NC_025457.alt2 NC_025457.ref  89.3701  508   34017 34524 22188 22695 454   54
NC_025457.alt2 NC_025457.ref  88.0240  501   18330 18830 11549 12049 441   60
```

| Column | Description |
| --- | --- |
| query | Identifier (name) of query sequence |
| reference | Identifier (name) of reference sequence |
| pident | Percent identity of local alignment |
| alnlen | Alignment length |
| qstart | Start of alignment in query |
| qend | End of alignment in query |
| rstart | Start of alignment in reference |
| rend | End of alignment in reference |
| nt_match | Number of matched (identical) nucleotides  |
| nt_mismatch | Number of mismatching nucleotides |


## Further clustering

The LZ-ANI output files, [ani.tsv](./example/output/ani.tsv) and [ani.ids.tsv](./example/output.ani.ids.tsv), can be used as input for clustering with [Clusty](https://github.com/refresh-bio/clusty). Clustering can use one of similarity measures (e.g., `tani`, `ani`), with the user specifying the minimum similarity threshold for connecting genomes.

For example, to cluster genomes with tANI ≥ 0.95, use the following command:

```bash
clusty --objects-file example/output/ani.ids.tsv --algo complete --distance-col tani --similarity --numeric-ids --min tani 0.95 example/output/ani.tsv clusters.txt
```

Clusty can also apply additional thresholds for various similarity measures. If a genome pair does not meet these thresholds, it is excluded from clustering.

```bash
# Cluster genomes based on ANI, connecting them only if ANI ≥ 95% and coverage ≥ 85%.
clusty --objects-file example/output/ani.ids.tsv --algo complete --distance-col ani --similarity --numeric-ids --min ani 0.95 --min qcov 0.85 example/output/ani.tsv clusters.txt
```

## Cite

Zielezinski A, Gudyś A, Barylski J, Siminski K, Rozwalak P, Dutilh BE, Deorowicz S. *Ultrafast and accurate sequence alignment and clustering of viral genomes*. bioRxiv [[doi](https://google.pl)][[pubmed](https://google.pl)].
