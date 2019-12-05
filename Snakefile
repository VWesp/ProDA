configfile: "config.yaml"

rule finish:
    input:
        "retained_finished.txt",
        "discarded_finished.txt"

include: "snakefiles/BLAST.snakefile"

include: "snakefiles/QUERY_CONTIG_MATCHER.snakefile"

include: "snakefiles/EXONERATE.snakefile"

include: "snakefiles/SPALN.snakefile"

include: "snakefiles/CD_HIT.snakefile"

include: "snakefiles/ALIGNMENT.snakefile"

include: "snakefiles/JOIN_RESULTS.snakefile"

include: "snakefiles/VISUALIZE_RETAINED_RESULTS.snakefile"

include: "snakefiles/VISUALIZE_DISCARDED_RESULTS.snakefile"
