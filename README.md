# Immuno

Immuno is a peptide-scoring pipeline designed to identify highly immunogenic peptide candidates from genomic or proteomic sequences. The tool computes an integrated score based on:

- **Immunogenicity** (allele-independent, fast to compute) [immunogenicity.py](https://github.com/comp-bio/immuno/blob/main/app/core/immunogenicity.py)
- **MHC-I binding affinity** (allele-specific, slower to compute) [iedb.org](https://downloads.iedb.org/tools/mhci/)

By filtering candidates based on immunogenicity first, Immuno significantly reduces the number of affinity evaluations required.

## Key Features

- Generates all possible peptides from a FASTA sequence
- Computes immunogenicity for each peptide
- Applies adaptive thresholding to skip affinity evaluation for weak candidates
- Integrates with IEDB MHC-I tools for affinity prediction
- Produces an integrated score and returns top candidates

## Algorithm Overview

The pipeline evaluates all possible peptides extracted from the input FASTA file.
To reduce the number of affinity predictions, Immuno applies the following strategy.

### 1. Immunogenicity computation and sorting

For every peptide, the immunogenicity (IM) score is computed first.
All peptides are then sorted in descending order of immunogenicity.

### 2. Iterative affinity evaluation

Affinity (AF) is computed only when needed, starting from peptides with the highest immunogenicity.

### 3. Early stopping

As peptides are processed in decreasing immunogenicity order, the algorithm keeps collecting peptides whose final integrated score satisfies the user-defined ranking criteria. If the required number of top peptides (defined by -n) is reached, the algorithm stops immediately.

### 4. Full evaluation fallback

If the algorithm reaches the end of the sorted list without obtaining the required number of top-scoring peptides:
either the input sequence does not contain enough suitable candidates or the requested number of results (-n) is too large.
In this case, all affinities will be computed, and the available candidates are returned.

## Installation

### IEDB MHC-I Tools

**Prerequisites**

- Linux 64-bit
- Python ≥ 3.6
- tcsh
- gawk

```bash
sudo apt-get install tcsh
sudo apt-get install gawk
```

**Install IEDB MHC-I**

```bash
wget https://downloads.iedb.org/tools/mhci/LATEST/IEDB_MHC_I-3.1.6.tar.gz
tar -zxvf IEDB_MHC_I-3.1.6.tar.gz
cd mhc_i
./configure
```

## Usage

The pipeline computes immunogenicity and affinity for all peptides generated from the input FASTA file and produces a ranked list of candidates.

### Basic command

```
run.py -src <fasta> [-freq <file> OR -allele <Name>] [options]
```

## Options

```
-src       <fasta-file>
-freq      <tsv-file with allele and frequency>
-allele    <allele name, e.g. "HLA-A*23:01">
-organism  [mouse | human(default)]
-model     [netmhcpan_el | netmhcpan_ba | consensus(default)]
-dir       <output directory>
-n         number of top results to return (10–1000, default 20)
-t         thread count
```

## Examples

```bash
run.py -src ./data/n5.fa -freq ./data/freq.tsv -model netmhcpan_ba
run.py -src ./data/na.fa -allele "HLA-A*30:01" -dir ./output/NA-res
run.py -src ./data/na.fa -freq ./data/freq.tsv
```

## Output

The tool generates:

- Ranked peptide list
- Integrated immunogenicity/affinity score
- Affinity predictions for selected alleles
- Summary statistics
