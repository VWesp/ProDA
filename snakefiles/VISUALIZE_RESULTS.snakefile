import os

rule visualize_results:
    input:
        "results/retained/proda.tsv",
        "results/discarded/proda.tsv"
    output:
        "results/retained/alignments/{subject}/{query}.png",
        temp("results/retained/alignments/{subject}/{query}_hit.faa"),
        temp("results/retained/alignments/{subject}/{query}_query.faa"),
        temp("results/retained/alignments/{subject}/{query}.clustal"),
        "results/discarded/alignments/{subject}/{query}.png",
        temp("results/discarded/alignments/{subject}/{query}_hit.faa"),
        temp("results/discarded/alignments/{subject}/{query}_query.faa"),
        temp("results/discarded/alignments/{subject}/{query}.clustal")
    params:
        props=config["properties"]
    log:
        "log/results/retained/alignments/{subject}/{query}.log",
        "log/results/discarded/alignments/{subject}/{query}.log"
    run:
        try:
            with open(input[0], "r") as proda_reader:
                next(proda_reader)
                content = proda_reader.readlines()
                for line in content:
                    subject = line.split("\t")[0]
                    query = line.split("\t")[1]
                    hit_sequence = line.split("\t")[4]
                    query_sequence = line.split("\t")[5]
                    with open(output[1], "w") as hit_writer:
                        hit_writer.write(">" + subject + "\n" + hit_sequence)

                    with open(output[2], "w") as query_writer:
                        query_writer.write(">" + query + "\n" + query_sequence)

                    os.system("(stretcher -asequence " + output[1] + " -sprotein1 -bsequence " +
                              output[2] + " -auto -aformat clustal -stdout > " + output[3] + ") 2> " + log[0])

                    os.system("(jalview -nodisplay -props " + params[0] + " -colour clustal -open " +
                              output[3] + " -png " + output[0] + ") 2> " + log[0])
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(ex))

        try:
            with open(input[1], "r") as proda_reader:
                next(proda_reader)
                content = proda_reader.readlines()
                for line in content:
                    subject = line.split("\t")[0]
                    query = line.split("\t")[1]
                    hit_sequence = line.split("\t")[4]
                    query_sequence = line.split("\t")[5]
                    with open(output[5], "w") as hit_writer:
                        hit_writer.write(">" + subject + "\n" + hit_sequence)

                    with open(output[6], "w") as query_writer:
                        query_writer.write(">" + query + "\n" + query_sequence)

                    os.system("(stretcher -asequence " + output[5] + " -sprotein1 -bsequence " +
                              output[6] + " -auto -aformat clustal -stdout > " + output[7] + ") 2> " + log[1])

                    os.system("(jalview -nodisplay -props " + params[0] + " -colour clustal -open " +
                              output[7] + " -png " + output[4] + ") 2> " + log[1])
        except Exception as ex:
            with open(log[1], "w") as log_writer:
                log_writer.write(str(ex))
