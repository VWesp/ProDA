import os

configfile: "config.yaml"

rule all:
    input:
        "temp/end.txt"

include: "snakefiles/BLAST.snakefile"

rule matcher:
    input:
        expand("data/subjects/{subject}.fna", subject=config["subjects"]),
        expand("blast_results/{subject}/{query}.hit", subject=config["subjects"], query=config["queries"])
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

rule test:
    input:
        expand("cd_hit/exonerate/{subject}/{query}.faa", subject=config["subjects"], query=config["queries"]),
        expand("cd_hit/spaln/{subject}/{query}.faa", subject=config["subjects"], query=config["queries"])
    output:
        "temp/end.txt"
    shell:
        "touch {output}"
