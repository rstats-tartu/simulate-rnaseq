from snakemake.utils import min_version
from Bio.SeqIO.FastaIO import SimpleFastaParser
from os.path import dirname
import numpy as np
import pandas as pd

min_version("5.1.2")


configfile: "config.yaml"


include: "rules/common.smk"


WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers"

# Number of samples
N = config["replicates"]

# Number of groups
N_GROUPS = config["n_groups"]

CONTRASTS = config["contrasts"]
SAMPLES = [f'sample_{"{:0>2}".format(i)}' for i in list(range(1, N_GROUPS * N + 1, 1))]


rule all:
    input:
        expand(
            [
                "output/diffexp/{contrast}.diffexp.tsv",
                "output/diffexp/{contrast}.ma-plot.svg",
            ],
            contrast=CONTRASTS,
        ),


report: "report/workflow.rst"


localrules:
    all,
    transcripts,


# Downlads human protein transcripts
rule transcripts:
    output:
        "data/transcripts.fa",
    log:
        "logs/transcripts.log",
    shell:
        """
        (mkdir -p data && wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.pc_transcripts.fa.gz | gunzip -c > {output[0]}) 2> {log}
        """


N_TX = config["n_transcripts"]["n"]
SEED_TX = config["n_transcripts"]["seed"]


rule sample_transcripts:
    """
    Sample desired number of transcripts from whole transcriptome
    """
    input:
        input=rules.transcripts.output[0],
    output:
        out="output/transcripts_sample.fasta",
    log:
        "logs/sample_transcripts.log",
    params:
        extra=f"samplereadstarget={N_TX} sampleseed={SEED_TX}",
    resources:
        runtime=120,
        mem_mb=4000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.6/bbtools/reformat"


rule simulate_experiment:
    input:
        fasta=rules.sample_transcripts.output.out,
    output:
        simulated_reads=expand(
            "output/simulated_reads/{sample}_{pair}.fasta", sample=SAMPLES, pair=[1, 2]
        ),
        samples="output/simulated_reads/sim_rep_info.txt",
    log:
        "logs/simulate_reads.log",
    params:
        reads_per_transcript=lambda wildcards, input: reads_per_transcript(
            input.fasta, 20
        ),
        fold_changes=fold_changes(
            [0.5, 1, 2],
            prob=[0.05, 0.9, 0.05],
            n_groups=N_GROUPS,
            n_transcripts=N_TX,
            seed=config["fold_changes"]["seed"],
        ),
        num_reps=np.repeat(N, N_GROUPS, axis=0),
        outdir=lambda wildcards, output: dirname(output.simulated_reads[0]),
        seed=config["simulate_experiment"]["seed"],
    resources:
        mem_mb=4000,
        runtime=120,
    wrapper:
        "file:wrappers/polyester"


rule shuffle:
    input:
        r1="output/simulated_reads/{sample}_1.fasta",
        r2="output/simulated_reads/{sample}_2.fasta",
    output:
        sh_r1="output/simulated_reads/shuffled_{sample}_1.fasta",
        sh_r2="output/simulated_reads/shuffled_{sample}_2.fasta",
    log:
        "logs/shuffle/{sample}.log",
    conda:
        "envs/bbmap.yaml"
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        """
        shuffle.sh in={input.r1} in2={input.r2} out={output.sh_r1} out2={output.sh_r2} -Xmx{resources.mem_mb}m 2> {log}
        """


rule salmon_index:
    input:
        rules.sample_transcripts.output.out,
    output:
        directory("output/salmon/index"),
    log:
        "logs/salmon_index.log",
    params:
        # optional parameters
        extra="",
    threads: 8
    resources:
        mem_mb=8000,
        runtime=120,
    wrapper:
        "v0.69.0/bio/salmon/index"


rule salmon_quant_reads:
    input:
        # If you have multiple fastq files for a single sample (e.g. technical replicates)
        # use a list for r1 and r2.
        r1="output/simulated_reads/shuffled_{sample}_1.fasta",
        r2="output/simulated_reads/shuffled_{sample}_2.fasta",
        index="output/salmon/index",
    output:
        quant="output/salmon/{sample}/quant.sf",
        lib="output/salmon/{sample}/lib_format_counts.json",
    log:
        "logs/salmon/{sample}.log",
    params:
        # optional parameters
        libtype="A",
        #zip_ext = bz2 # req'd for bz2 files ('bz2'); optional for gz files('gz')
        extra="",
    threads: 8
    resources:
        mem_mb=8000,
        runtime=120,
    wrapper:
        "v0.69.0/bio/salmon/quant"


rule count_matrix:
    input:
        expand(
            "output/salmon/{sample}/quant.sf", sample=SAMPLES,
        ),
    output:
        "output/counts.tsv",
    log:
        "logs/count-matrix.log",
    params:
        samples=SAMPLES,
    conda:
        "envs/pandas.yaml"
    resources:
        mem_mb=4000,
        runtime=120,
    script:
        "scripts/count-matrix.py"


rule deseq2_init:
    input:
        counts="output/counts.tsv",
        samples="output/simulated_reads/sim_rep_info.txt",
    output:
        "output/deseq2.rds",
    conda:
        "envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log",
    threads: get_deseq2_threads()
    resources:
        mem_mb=4000,
        runtime=120,
    script:
        "scripts/deseq2-init.R"


rule deseq2:
    input:
        "output/deseq2.rds",
    output:
        table=report("output/diffexp/{contrast}.diffexp.tsv", "report/diffexp.rst"),
        ma_plot=report("output/diffexp/{contrast}.ma-plot.svg", "report/ma.rst"),
    params:
        contrast=get_contrast,
    conda:
        "envs/deseq2.yaml"
    log:
        "logs/{contrast}.diffexp.log",
    threads: get_deseq2_threads
    resources:
        mem_mb=4000,
        runtime=120,
    script:
        "scripts/deseq2.R"
