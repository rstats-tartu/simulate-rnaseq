
# Number of samples
N = 10
N_GROUPS = 2
SAMPLES = [ f'sample_{"{:0>2}".format(i)}' for i in list(range(1, N_GROUPS * N + 1, 1))]

rule all:
    input: "output/counts.tsv", expand("output/salmon/{sample}/quant.sf", sample = SAMPLES)


rule transcripts:
     output:
         "output/gencode.v36.pc_transcripts.fa",
     log:
         "logs/transcripts.log"
     shell:
         """
         (wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.pc_transcripts.fa.gz \
         && mkdir -p output \
         && mv gencode.v36.pc_transcripts.fa.gz output \
         && gunzip gencode.v36.pc_transcripts.fa.gz) 2> {log}
         """

rule simulate_reads:
    input:
        transcripts="output/gencode.v36.pc_transcripts.fa",
    output:
        transcripts_sample="output/transcripts_sample.fa",
        simulated_reads=expand("output/simulated_reads/{sample}_{pair}.fasta", sample = SAMPLES, pair = [1, 2]),
    log:
        "logs/simulate_reads.log",
    params:
        n_transcripts = 20000,
        seed = [11, 12, 13],
        fold_change_values=[0.5, 1, 2],
        n_groups = N_GROUPS,
        probs = [0.0025, 0.995, 0.0025],
        coverage = 20,
    conda:
        "envs/polyester.yaml"
    resources:
       mem_mb = 4000,
       runtime = 120,
    script:
       "scripts/simulate_reads.R"


rule shuffle:
    input:
        r1 = "output/simulated_reads/{sample}_1.fasta",
        r2 = "output/simulated_reads/{sample}_2.fasta",
    output:
        sh_r1 = "output/simulated_reads/shuffled_{sample}_1.fasta",
        sh_r2 = "output/simulated_reads/shuffled_{sample}_2.fasta",
    log:
        "logs/shuffle/{sample}.log"
    conda:
       "envs/bbmap.yaml"
    resources:
       mem_mb = 4000,
       runtime = 120,
    shell:
       """
       shuffle.sh in={input.r1} in2={input.r2} out={output.sh_r1} out2={output.sh_r2} -Xmx{resources.mem_mb}m 2> {log}
       """


rule salmon_index:
    input:
        "output/fasta_sample.fa"
    output:
        directory("output/salmon/index")
    log:
        "logs/salmon_index.log"
    params:
        # optional parameters
        extra=""
    threads: 8
    resources:
        mem_mb=8000,
        runtime=120
    wrapper:
        "v0.69.0/bio/salmon/index"

rule salmon_quant_reads:
    input:
        # If you have multiple fastq files for a single sample (e.g. technical replicates)
        # use a list for r1 and r2.
        r1 = "output/simulated_reads/shuffled_{sample}_1.fasta",
        r2 = "output/simulated_reads/shuffled_{sample}_2.fasta",
        index = "output/salmon/index"
    output:
        quant = "output/salmon/{sample}/quant.sf",
        lib = "output/salmon/{sample}/lib_format_counts.json"
    log:
        'logs/salmon/{sample}.log'
    params:
        # optional parameters
        libtype ="A",
        #zip_ext = bz2 # req'd for bz2 files ('bz2'); optional for gz files('gz')
        extra=""
    threads: 8
    resources:
        mem_mb=8000,
        runtime=120
    wrapper:
        "v0.69.0/bio/salmon/quant"


rule count_matrix:
    input:
        expand(
            "output/salmon/{sample}/quant.sf",
            sample=SAMPLES,
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
        runtime=120
    script:
        "scripts/count-matrix.py"

