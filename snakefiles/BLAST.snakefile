rule makeblastdb:
    input:
        subjects = lambda wildcards: config["subjects"][wildcards.subject]
    output:
        "blast_dbs/{subject}/temp.db"
    log:
        "log/blast_dbs/{subject}.log"
    shell:
        "(makeblastdb -in {input.subjects} -dbtype nucl "
        "-out blast_dbs/{wildcards.subject}/{wildcards.subject}) 2> {log} && "
        "touch {output}"

rule blast:
    input:
        temp_database = "blast_dbs/{subject}/temp.db",
        queries = lambda wildcards: config["queries"][wildcards.query]
    output:
        "blast_results/{subject}/{query}.hit"
    params:
        evalue=config["blast_evalue"]
    log:
        "log/blast_results/{subject}/{query}.log"
    shell:
        "(tblastn -query {input.queries} -db blast_dbs/{wildcards.subject}/{wildcards.subject} "
        "-outfmt '6 qseqid sseqid sstart send' "
        "-evalue {params.evalue} > {output}) 2> {log}"
