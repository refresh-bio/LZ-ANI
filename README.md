# LZ-ANI
LZ-ANI is a fast and memory-efficient tool for determining average nucleotide identity (ANI) in large sets of genomic sequences. The tool uses Lempel-Ziv parsing, achieving high sensitivity in discerning matched and mismatched nucleotides, thereby enabling accurate determination of ANI. Its efficiency stems from a simplified indel handling model, making LZ-ANI magnitudes faster than alignment-based tools (e.g., BLASTn, MegaBLAST) while maintaining accuracy comparable to the most sensitive BLASTn searches. 

LZ-ANI is a key component of [Vclust](https://github.com/refresh-bio/vclust), a tool for comparing and clustering virus genomes. Although optimized for virus genomes, LZ-ANI's parameters can be customized for longer genomes, such as those of bacteria and archaea.

#### Sequence similarity measures

LZ-ANI offers six similarity measures between two genomic sequences:

- **ANI**: The number of identical bases across local alignments divided by the total length of the alignments.
- **Global ANI (gANI)**: The number of identical bases across local alignments divided by the length of the query/target genome.
- **Total ANI (tANI)**: The number of identical bases between query-target and target-query genomes divided by the sum length of both genomes.
- **Coverage (alignment fraction)**: The proportion of the query sequence aligned with the target sequence.
- **Number of local alignments**: The count of individual alignments found between the sequences.
- **Ratio between query and target genome lengths**: A measure comparing the lengths of the two genomes.


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

For brevity, only the first 15 lines of output are shown:

```
id1 id2 tani  gani  ani cov len_ratio
NC_025457.alt2  NC_005091.alt2  0.013765  0.011564  0.577882  0.020011  1.007347
NC_005091.alt2  NC_025457.alt2  0.013765  0.015982  0.575792  0.027757  0.992706
NC_025457.alt2  NC_005091.alt1  0.014603  0.013995  0.565491  0.024749  1.116770
NC_005091.alt1  NC_025457.alt2  0.014603  0.015282  0.555345  0.027517  0.895440
NC_025457.alt2  NC_005091.ref   0.014644  0.012671  0.576596  0.021975  1.116770
NC_005091.ref   NC_025457.alt2  0.014644  0.016848  0.569077  0.029606  0.895440
NC_025457.alt2  NC_002486.alt   0.022687  0.018328  0.604938  0.030297  1.405995
NC_002486.alt   NC_025457.alt2  0.022687  0.028815  0.594216  0.048492  0.711240
NC_025457.alt2  NC_002486.ref   0.020692  0.017268  0.604474  0.028567  1.405995
NC_002486.ref   NC_025457.alt2  0.020692  0.025506  0.609424  0.041853  0.711240
NC_025457.alt2  NC_025457.ref   0.752589  0.658220  0.910059  0.723272  1.504290
NC_025457.ref   NC_025457.alt2  0.752589  0.894547  0.915166  0.977470  0.664765
NC_025457.alt2  NC_025457.alt1  0.595191  0.502322  0.895679  0.560829  1.562460
NC_025457.alt1  NC_025457.alt2  0.595191  0.740296  0.909148  0.814275  0.640016
NC_025457.alt2  NC_010807.alt2  0.027875  0.022115  0.570567  0.038760  1.582148
```

### Output format

The `--out-format` provides three output views: `standard`, `lite`, and `complete`.

| Field | Standard | Lite | Complete | Description |
| --- | :---: |:---: | :---: | --- |
| idx1 | + | + | +  | index of sequence 1 |
| idx2 | + | + | +  | index of sequence 2 |
| id1 | + | - | +  | identifier (name) of sequence 1 |
| id2 | + | - | +  | identifier (name) of sequence 2 |
| tani | + | + | +  | total ANI [0-1] |
| gani | + | + | +  | global ANI [0-1] |
| ani | + | + | +  | ANI [0-1] |
| cov | + | + | +  | Coverage (alignment fraction) [0-1] |
| num_alns | + | + | +  | Number of alignments |
| len_ratio | + | + | +  | Length ratio between sequence 1 and sequence 2 |
| len1 | - | - | +  | Length of sequence 1 |
| len2 | - | - | +  | Length of sequence 2|
| nt_match | - | - | +  | Number of matching nucleotides across alignments |
| nt_mismatch | - | - | +  | Number of mismatching nucleotides across alignments |


In addition, the `--out-format` option permits formatting arbitrary fields from the LZ-ANI tab-separated-value (TSV) format: 

```bash
./lz-ani all2all --in-fasta example/multifasta.fna --out ani.tsv --out-format id1,id2,ani,cov
```

```
id1 id2 ani cov
NC_025457.alt2  NC_005091.alt2  0.577882  0.020011
NC_005091.alt2  NC_025457.alt2  0.575792  0.027757
NC_025457.alt2  NC_005091.alt1  0.565491  0.024749
NC_005091.alt1  NC_025457.alt2  0.555345  0.027517
NC_025457.alt2  NC_005091.ref   0.576596  0.021975
NC_005091.ref   NC_025457.alt2  0.569077  0.029606
NC_025457.alt2  NC_002486.alt   0.604938  0.030297
NC_002486.alt   NC_025457.alt2  0.594216  0.048492
NC_025457.alt2  NC_002486.ref   0.604474  0.028567
NC_002486.ref   NC_025457.alt2  0.609424  0.041853
NC_025457.alt2  NC_025457.ref   0.910059  0.723272
NC_025457.ref   NC_025457.alt2  0.915166  0.977470
NC_025457.alt2  NC_025457.alt1  0.895679  0.560829
NC_025457.alt1  NC_025457.alt2  0.909148  0.814275
NC_025457.alt2  NC_010807.alt2  0.570567  0.038760
...
```


### Output filtering

The `--out-filter` option allows you to filter the output by setting minimum similarity thresholds, enabling you to report only those genome pairs that meet the specified criteria, thus significantly reducing the output TSV file size. For example, the following command outputs only genome pairs with ANI ≥ 0.95 and coverage ≥ 0.85:

```bash
./lz-ani all2all --in-fasta example/multifasta.fna --out ani.tsv --out-filter ani 0.95 --out-filter cov 0.85
```

```
id1 id2 tani  gani  ani cov len_ratio
NC_005091.alt2  NC_005091.ref   0.966298  0.967989  0.968125  0.999859  1.108624
NC_005091.ref   NC_005091.alt2  0.966298  0.964424  0.964911  0.999495  0.902019
NC_005091.alt1  NC_005091.ref   0.970072  0.970151  0.971368  0.998747  1.000000
NC_005091.ref   NC_005091.alt1  0.970072  0.969994  0.971245  0.998712  1.000000
NC_002486.alt   NC_002486.ref   1.000000  1.000000  1.000000  1.000000  1.000000
NC_002486.ref   NC_002486.alt   1.000000  1.000000  1.000000  1.000000  1.000000
NC_025457.alt1  NC_025457.ref   0.809496  0.845785  0.985613  0.858131  0.962770
NC_010807.alt2  NC_010807.alt3  0.972839  0.985279  0.987642  0.997608  1.016645
NC_010807.alt3  NC_010807.alt2  0.972839  0.960192  0.986657  0.973177  0.983627
NC_010807.alt2  NC_010807.alt1  0.987250  0.987449  0.987547  0.999901  1.044828
NC_010807.alt1  NC_010807.alt2  0.987250  0.987041  0.987117  0.999923  0.957095
NC_010807.alt2  NC_010807.ref   0.989807  0.990063  0.990063  1.000000  1.044828
NC_010807.ref   NC_010807.alt2  0.989807  0.989540  0.989617  0.999923  0.957095
NC_010807.alt3  NC_010807.alt1  0.979963  0.967035  0.994304  0.972575  1.027721
NC_010807.alt1  NC_010807.alt3  0.979963  0.993250  0.994557  0.998686  0.973026
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

## Further clustering

The LZ-ANI output files, [ani.tsv](./example/output/ani.tsv) and [ani.ids.tsv](./example/output.ani.ids.tsv), can be used as input for clustering with [Clusty](https://github.com/refresh-bio/clusty). Clustering can use one of similarity measures (e.g., `tani`, `ani`), with the user specifying the minimum similarity threshold for connecting genomes.

For example, to cluster genomes with tANI ≥ 0.95, use the following command:

```bash
clusty --objects-file example/output/ani.ids.tsv --algo complete --distance-col tani --similarity --numeric-ids --min tani 0.95 example/output/ani.tsv clusters.txt
```

Clusty can also apply additional thresholds for various similarity measures. If a genome pair does not meet these thresholds, it is excluded from clustering.

```bash
# Cluster genomes based on ANI, connecting them only if ANI ≥ 95% and coverage ≥ 85%.
clusty --objects-file example/output/ani.ids.tsv --algo complete --distance-col ani --similarity --numeric-ids --min ani 0.95 --min cov 0.85 example/output/ani.tsv clusters.txt
```

## Cite

Zielezinski A, Gudyś A, Barylski J, Siminski K, Rozwalak P, Dutilh BE, Deorowicz S. *Ultrafast and accurate sequence alignment and clustering of viral genomes*. bioRxiv [[doi](https://google.pl)][[pubmed](https://google.pl)].
