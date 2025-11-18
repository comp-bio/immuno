# Immuno

## Installation

### MHC_I

Prerequisites:

1. Linux 64-bit environment
2. Python 3.6 or higher
3. tcsh http://www.tcsh.org/Welcome (`sudo apt-get install tcsh`)
4. gawk http://www.gnu.org/software/gawk/ (`sudo apt-get install gawk`)

> https://downloads.iedb.org/tools/mhci/LATEST/

```bash
wget https://downloads.iedb.org/tools/mhci/LATEST/IEDB_MHC_I-3.1.6.tar.gz
tar -zxvf IEDB_MHC_I-3.1.6.tar.gz
cd mhc_i
./configure
```

### Immuno

The pipeline calculates a metric integrated from affinity and
immunogenicity for all possible peptides from the specified .fasta file

Usage:

```text
  run.py -src <fasta> [-freq <file> OR -allele [Name]] [options]
```

Options:

```text
  -src       <fasta-file>
  -freq      <tsv-file with two columns: allele name and frequency>
  -allele    [allele name, for example "HLA-A*23:01"]
  -organism  [organism: mouse, human (default)]
  -model     [model name: netmhcpan_el, netmhcpan_ba, consensus (default)]
  -dir       <results dir name>
  -n         [number] results count, 10–1000 (default 20)
  -dist      [number] Distance between peptides for clustering, 0–1 (default 0.15)
  -t         [number] thread count
```

Examples:

```bash
  run.py -src ./data/n5.fa -freq ./data/freq.tsv -model netmhcpan_ba
  run.py -src ./data/na.fa -allele "HLA-A*30:01" -dir ./output/NA-res
  run.py -src ./data/na.fa -freq ./data/freq.tsv
```
