def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(SAMPLES) < 100 or few_coeffs else 6


def get_contrast(wildcards):
    return CONTRASTS[wildcards.contrast]


def reads_per_transcript(input, coverage):
    seqlen = []
    with open(input, "rU") as fastafile:
        for name, seq in SimpleFastaParser(fastafile):
            seqlen.append(len(seq))
    return [round(i * coverage / 100) for i in seqlen]


def fold_changes(values, prob, n_groups, n_transcripts, seed=11):
    np.random.seed(seed)
    return np.random.choice(values, size=n_groups * n_transcripts, replace=True, p=prob)
