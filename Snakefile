configfile: "config.yaml"

rule all:
    input:
        "finished.txt"

include: "snakefiles/BLAST.snakefile"

include: "snakefiles/QUERY_CONTIG_MATCHER.snakefile"

include: "snakefiles/EXONERATE.snakefile"

include: "snakefiles/SPALN.snakefile"

include: "snakefiles/CD_HIT.snakefile"

include: "snakefiles/ALIGNMENT.snakefile"

include: "snakefiles/GET_BEST_HIT.snakefile"

include: "snakefiles/JOIN_RESULTS.snakefile"

include: "snakefiles/VISUALIZE_RESULTS.snakefile"

rule finisher:
    input:
        #expand("results/retained/alignments/{subject}/{query}.png", subject=config["subjects"], query=config["queries"]),
        #expand("results/discarded/alignments/{subject}/{query}.png", subject=config["subjects"], query=config["queries"])
        "results/retained/proda.tsv",
        "results/discarded/proda.tsv"
    output:
        temp("finished.txt")
    shell:
        "touch {output}"
