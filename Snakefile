configfile: "config.yaml"

rule finish:
    input:
        expand("temp/above/{subject}/{query}.txt", subject=config["subjects"], query=config["queries"]),
        expand("temp/below/{subject}/{query}.txt", subject=config["subjects"], query=config["queries"])

include: "snakefiles/BLAST.snakefile"

include: "snakefiles/QUERY_CONTIG_MATCHER.snakefile"

include: "snakefiles/EXONERATE.snakefile"

include: "snakefiles/EX_GFF_PARSER.snakefile"

include: "snakefiles/EX_STRETCHER.snakefile"

include: "snakefiles/SPALN.snakefile"

include: "snakefiles/SP_GFF_PARSER.snakefile"

include: "snakefiles/SP_STRETCHER.snakefile"

include: "snakefiles/JOIN_RESULTS.snakefile"

include: "snakefiles/THRESHOLD_FILTER.snakefile"

include: "snakefiles/VISUALIZE_RETAINED_RESULTS.snakefile"

include: "snakefiles/VISUALIZE_DISCARDED_RESULTS.snakefile"
