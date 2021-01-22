log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(polyester)
library(Biostrings)

tx <- snakemake@input[["transcripts"]]
tx_sample <- snakemake@output[["transcripts_sample"]]
num_reps <- snakemake@params[["n_reps"]]
n_groups <- snakemake@params[["n_groups"]]
tx_n <- snakemake@params[["n_transcripts"]]
fold_change_values <- snakemake@params[["fold_change_values"]]
prob <- snakemake@params[["fold_change_probs"]]

tx_count <- count_transcripts(tx)
fasta <- readDNAStringSet(tx)

set.seed(snakemake@params[["seed"]][1])
fasta_sample_idx <- sample(1:tx_count, size = tx_n)
fasta_sample <- fasta[fasta_sample_idx]
writeXStringSet(fasta_sample, tx_sample)

set.seed(snakemake@params[["seed"]][2])
fold_change_values <- sample(
  fold_change_values,
  size = n_groups * tx_n,
  prob = prob, 
  replace = TRUE
)

fold_changes <- matrix(fold_change_values, nrow = tx_n)
readspertx <- round(snakemake@params[["coverage"]] * width(fasta_sample) / 100)

# remove quotes from transcript IDs:
tx_names = gsub("'", "", names(readDNAStringSet(tx_sample)))

simulate_experiment(
  tx_sample, 
  reads_per_transcript = readspertx,
  fold_changes = fold_changes, 
  num_reps = rep(num_reps, n_groups),
  outdir = dirname(snakemake@output[["simulated_reads"]])[1],
  transcriptid = tx_names,
  seed = snakemake@params[["seed"]][3]
)

