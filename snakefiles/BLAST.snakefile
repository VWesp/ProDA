rule Build_DNA_Databases:
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

rule Find_Candidates:
    input:
        "blast_dbs/{subject}/temp.db",
        "data/queries/{query}.faa"
    output:
        "blast_results/{subject}/{query}.hit"
    params:
        evalue=config["blast_evalue"]
    threads: config["threads"]
    log:
        "log/blast_results/{subject}/{query}.log"
    shell:
        "(tblastn -query {input[1]} -db blast_dbs/{wildcards.subject}/{wildcards.subject} "
        "-outfmt '6 qseqid sseqid sstart send' "
        "-evalue {params.evalue} -num_threads {threads} > {output}) 2> {log}"
