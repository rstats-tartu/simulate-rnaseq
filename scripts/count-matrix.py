
import sys
sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

counts = [pd.read_table(f, index_col=0, usecols=[0, 4], header=None, skiprows=4) for f in snakemake.input]
for t, sample in zip(counts, snakemake.params.samples):
    t.columns = [sample]
matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
matrix = matrix.groupby(matrix.columns, axis=1).sum()

matrix.to_csv(snakemake.output[0], sep="\t")

