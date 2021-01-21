log <- file(snakemake@log[[1]], open="wt")
sink(log, type="message")

library(polyester)
library(Biostrings)

tx_count <- count_transcripts(snakemake@input[["transcripts"]])
fasta <- readDNAStringSet(snakemake@input[["transcripts"]])
set.seed(snakemake@params[["seed"]][1])
tx_n <- snakemake@params[["tx_n"]]
fasta_sample_idx <- sample(1:tx_count, size = tx_n)
fasta_sample <- fasta[fasta_sample_idx]
writeXStringSet(fasta_sample, snakemake@output[["transcripts_sample"]])

set.seed(snakemake@params[["seed"]][2])
fold_change_values <- sample(
  snakemake@params[["fold_change_values"]], 
  size = snakemake@params[["n_groups"]] * tx_n,
  prob = snakemake@params[["fold_change_probs"]], 
  replace = TRUE
)

fold_changes <- matrix(fold_change_values, nrow = tx_n)
readspertx <- round(snakemake@params[["coverage"]] * width(fasta_sample) / 100)

# remove quotes from transcript IDs:
tx_names = gsub("'", "", names(readDNAStringSet(snakemake@output[["transcripts_sample"]])))

simulate_experiment(
  snakemake@output[["transcripts_sample"]], 
  reads_per_transcript = readspertx,
  fold_changes = fold_changes, 
  outdir = dirname(snakemake@output[["simulated_reads"]],
  transcriptid = tx_names,
  seed = snakemake@params[["seed"]][3]
)

