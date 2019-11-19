configfile: "config.yaml"

rule all:
    input:
        "temp/end.txt"

include: "snakefiles/BLAST.snakefile"

rule matcher:
    input:
        expand("data/subjects/{subject}.fasta", subject=config["subjects"]),
        expand("blast_results/{subject}/{query}.hit", subject=config["subjects"], query=config["queries"])
    output:
        "matches/{subject}/{query}.fasta"
    log:
        "log/matches/{subject}/{query}.log"
    params:
        left=config["left_addendum"],
        right=config["right_addendum"]
    script:
        "scripts/contig_query_matcher.py"

include: "snakefiles/EXONERATE.snakefile"

rule test:
    input:
        expand("exonerate_results/{subject}/{query}.fasta", subject=config["subjects"], query=config["queries"])
    output:
        "temp/end.txt"
    shell:
        "touch {output}"
