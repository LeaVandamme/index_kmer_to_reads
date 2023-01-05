import os
import itertools

# Path

path_log = "log"
path_reads = os.path.join(path_log,"reads")
path_sequences = os.path.join(path_log,"sequences")
path_memory_results = os.path.join(path_log,"memory_res")
path_reads_results = os.path.join(path_log,"reads_res")

reads_tests = ["simpleReads", "repetitions"]
reads_phage = ["phage10", "phage100", "phage1000"]
sequences_tests = ["unique", "repeted"]
sequences_phage = ["testPhage15"]
k_tests = [8]
k_phage = [15]
b = [26, 32]

# Les différentes exps

rule all:
    input:
        expand(os.path.join(path_memory_results,"{file}-{sequence}-{k}-{b}-memory-file.txt"),k=k_tests,b=b,file=reads_tests,sequence=sequences_tests),
        expand(os.path.join(path_memory_results,"{file}-{sequence}-{k}-{b}-memory-file.txt"),k=k_phage,b=b,file=reads_phage,sequence=sequences_phage),

# Les autres regles

rule creation_query:
    input:
        reads = os.path.join(path_reads, "{file}.fa"),
        sequences = os.path.join(path_sequences, "{sequence}.fa"),
    output:
        memory = os.path.join(path_memory_results,"{file}-{sequence}-{k}-{b}-memory-file.txt"),
        reads = os.path.join(path_reads_results,"{file}-{sequence}-{k}-{b}.fa"),

    shell:
        "make -j && ./index {input.reads} {wildcards.k} {wildcards.b} {output.reads} {output.memory} -s {input.sequences}"

