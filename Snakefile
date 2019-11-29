configfile: "config.yaml"

rule finish:
    input:
        "results/retained/proda.tsv",
        "results/discarded/proda.tsv",
        "retained_finished.txt",
        "discarded_finished.txt"

include: "snakefiles/BLAST.snakefile"

include: "snakefiles/QUERY_CONTIG_MATCHER.snakefile"

include: "snakefiles/EXONERATE.snakefile"

include: "snakefiles/SPALN.snakefile"

include: "snakefiles/CD_HIT.snakefile"

include: "snakefiles/ALIGNMENT.snakefile"

include: "snakefiles/GET_BEST_HIT.snakefile"

include: "snakefiles/JOIN_RESULTS.snakefile"

include: "snakefiles/VISUALIZE_RETAINED_RESULTS.snakefile"

include: "snakefiles/VISUALIZE_DISCARDED_RESULTS.snakefile"

include: "snakefiles/RETAINED_MAPPER.snakefile"

include: "snakefiles/DISCARDED_MAPPER.snakefile"

rule Join_Mapping_Results:
    input:
        expand("results/retained/proda_{subject}.tsv", subject=config["subjects"]),
        expand("results/discarded/proda_{subject}.tsv", subject=config["subjects"])
    output:
        "results/retained/proda.tsv",
        "results/discarded/proda.tsv"
    log:
        "log/results/retained/proda.log",
        "log/results/discarded/proda.log"
    run:
        retained_results = []
        discarded_results = []
        first_file = True
        for tsv in input[0]:
            if(os.stat(tsv).st_size != 0):
                with open(tsv, "r") as map_reader:
                    if(first_file):
                        retained_results.append(map_reader.read().strip())
                        first_file = False
                    else:
                        retained_results.append("\n".join(map_reader.read().split("\n")[1:]).strip())

        for tsv in input[0]:
            if(os.stat(tsv).st_size != 0):
                with open(tsv, "r") as map_reader:
                    if(first_file):
                        discarded_results.append(map_reader.read().strip())
                        first_file = False
                    else:
                        discarded_results.append("\n".join(map_reader.read().split("\n")[1:]).strip())

        with open(output[0], "w") as retained_writer:
            retained_writer.write("\n".join(retained_results))

        with open(output[1], "w") as discarded_writer:
            discarded_writer.write("\n".join(discarded_results))
