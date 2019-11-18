configfile: "config.yaml"

rule all:
    input:
        "temp/end.txt"

include: "snakefiles/BLAST.snakefile"

rule matcher:
    input:
        expand("subjects/{subject}.fasta", subject=config["subjects"]),
        expand("blast_results/{subject}/{query}.hit", subject=config["subjects"], query=config["queries"])
    output:
        "temp/scripts.txt",
    params:
        left=config["left_addendum"],
        right=config["right_addendum"]
    script:
        "scripts/contig_query_matcher.py"

rule test:
    input:
        "temp/scripts.txt"
    output:
        "temp/end.txt"
    shell:
        "touch {output}"
