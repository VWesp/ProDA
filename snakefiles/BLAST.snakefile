rule makeblastdb:
    input:
        "data/subjects/{subject}.fna"
    output:
        "blast_dbs/{subject}/temp.db"
    log:
        "log/blast_dbs/{subject}.log"
    shell:
        "(makeblastdb -in {input} -dbtype nucl "
        "-out blast_dbs/{wildcards.subject}/{wildcards.subject}) 2> {log} && "
        "touch {output}"

rule blast:
    input:
        "blast_dbs/{subject}/temp.db",
        "data/queries/{query}.faa"
    output:
        "blast_results/{subject}/{query}.hit"
    params:
        evalue=config["blast_evalue"],
    log:
        "log/blast_results/{subject}/{query}.log"
    shell:
        "(tblastn -query {input[1]} -db blast_dbs/{wildcards.subject}/{wildcards.subject} "
        "-outfmt '6 qseqid sseqid sstart send' "
        "-evalue {params.evalue} > {output}) 2> {log}"
