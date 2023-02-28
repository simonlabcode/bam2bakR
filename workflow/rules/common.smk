SAMP_NAMES = list(config['samples'].keys())

CTL_NAMES = list(config['control_samples'])

nctl = len(CTL_NAMES)

def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]

def get_index_name():
    genome = config["genome_fasta"]
    index = str(genome) + ".fai"
    return index