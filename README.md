# simulate-rnaseq

Simple workflow to simulate RNA-seq reads and perform DE-analysis.


Workflow simulates RNA-seq reads using **Polyester** R package <https://github.com/alyssafrazee/polyester> and 
performs DE-analysis using **DESeq2** R package <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>.

Human protein coding transcript sequences are downloaded from <https://www.gencodegenes.org/human/>. 
Transcript sequences can be customised by directly editing Snakefile.

# Acknowledgements

The workflow borrowed bits and pieces from rna-seq-star-deseq2 <https://github.com/snakemake-workflows/rna-seq-star-deseq2>.
