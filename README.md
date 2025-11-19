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
- Optional peptide clustering to avoid redundant peptides

## Algorithm Overview

Immunogenicity (IM) is computationally inexpensive, while affinity (AF) requires calling external tools.  
To optimize performance, Immuno uses an adaptive immunogenicity threshold $$T$$.

### 1. Initial filtering

A threshold $$T$$ is chosen (may be estimated using a binomial model).  
All peptides with $$IM < T$$ are skipped at this stage.

### 2. Affinity evaluation

Affinity is computed **only** for peptides with $$IM \ge T$$

### 3. Score upper bound

For any peptide with $$IM < T$$, the maximum possible integrated score is:  
$$\delta = \frac{1}{4}\left(2 - \sqrt{2}(1 - T) + 2T\right)$$

### 4. Early stopping

If enough peptides are found with $$score > \delta$$ the algorithm terminates early, because peptides with $$IM < T$$ cannot exceed this score.

### 5. Threshold adjustment

If not enough candidates exceed $$\delta$$, the threshold $$T$$ is reduced, and the process repeats.  
This technique reduces the total number of affinity evaluations, accelerating large-scale peptide screening.

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
-dist      clustering distance between peptides (0–1, default 0.15)
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
- Optional clustering results
- Summary statistics
