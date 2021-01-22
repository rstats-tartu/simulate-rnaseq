# simulate-rnaseq

Simple workflow to simulate RNA-seq reads and perform DE-analysis.


Workflow simulates RNA-seq reads using **Polyester** R package <https://github.com/alyssafrazee/polyester>, quantifies expression with **Salmon** <https://combine-lab.github.io/salmon/> and 
performs DE-analysis using **DESeq2** R package <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>.

Human protein coding transcript sequences are downloaded from <https://www.gencodegenes.org/human/>. 
Transcript sequences can be customised by directly editing Snakefile.

# Install

To install clone this repository and install dependencies (assumes that [conda/miniconda](https://docs.conda.io/en/latest/miniconda.html) is already installed):

```bash
conda env create -f environment.yml
```

This command creates conda environment called "simulaternaseq".


# Running

Activate conda environment:
```bash
conda activate simulaternaseq
```

For test run replace "-j" flag with "-n":
```bash

snakemake -j --use-conda
```


To generate snakemake run report:
```bash
snakemake --report output/report.html
```


# Some outputs of interest

File `output/diffexp/group1-vs-group2.diffexp.tsv` contains differential expression results, which can be compared to ground truth in `output/simulated_reads/sim_tx_info.txt`.


# Acknowledgements

The workflow borrowed bits and pieces from rna-seq-star-deseq2 <https://github.com/snakemake-workflows/rna-seq-star-deseq2>.
