import os

configfile: "config.yaml"

rule all:
    input:
        "temp/end.txt"

include: "snakefiles/BLAST.snakefile"

rule matcher:
    input:
        "data/subjects/{subject}.fna",
        "blast_results/{subject}/{query}.hit"
    output:
        "matches/{subject}/{query}.fna"
    log:
        "log/matches/{subject}/{query}.log"
    params:
        left=config["left_addendum"],
        right=config["right_addendum"]
    script:
        "scripts/contig_query_matcher.py"

include: "snakefiles/EXONERATE.snakefile"

include: "snakefiles/SPALN.snakefile"

include: "snakefiles/CD_HIT.snakefile"

include : "snakefiles/ALIGNMENT.snakefile"

include : "snakefiles/GET_BEST_HIT.snakefile"

rule test:
    input:
        expand("best_hit/{subject}/{query}.faa", subject=config["subjects"], query=config["queries"])
    output:
        temp("temp/end.txt")
    shell:
        "touch {output}"
