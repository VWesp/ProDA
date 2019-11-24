import os

rule visualize_results:
    input:
        "results/retained/proda.tsv",
        "results/discarded/proda.tsv"
    output:
        temp("finished.txt")
    params:
        props=config["properties"]
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
                    hit_output =  "results/retained/alignments/" + subject + "/" + query + "_hit.faa"
                    with open(hit_output, "w") as hit_writer:
                        hit_writer.write(">" + subject + "\n" + hit_sequence)

                    query_output =  "results/retained/alignments/" + subject + "/" + query + "_query.faa"
                    with open(query_output, "w") as query_writer:
                        query_writer.write(">" + query + "\n" + query_sequence)

                    clustal_output =  "results/retained/alignments/" + subject + "/" + query + ".clustal"
                    log_output = "log/results/retained/alignments/" + subject + "/" + query + ".log",
                    os.system("(stretcher -asequence " + hit_output + " -sprotein1 -bsequence " +
                              query_output + " -auto -aformat clustal -stdout > " + clustal_output + ") 2> " +log_output)

                    png_output = "results/retained/alignments/" + subject + "/" + query + ".png"
                    os.system("(jalview -nodisplay -props " + params[0] + " -colour clustal -open " +
                              clustal_output + " -png " + png_output + ") 2> " + log_output)
        except Exception as ex:
            with open("log/results/retained/alignments.log", "w") as log_writer:
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
                    hit_output =  "results/discarded/alignments/" + subject + "/" + query + "_hit.faa"
                    with open(hit_output, "w") as hit_writer:
                        hit_writer.write(">" + subject + "\n" + hit_sequence)

                    query_output =  "results/discarded/alignments/" + subject + "/" + query + "_query.faa"
                    with open(query_output, "w") as query_writer:
                        query_writer.write(">" + query + "\n" + query_sequence)

                    clustal_output =  "results/discarded/alignments/" + subject + "/" + query + ".clustal"
                    log_output = "log/results/discarded/alignments/" + subject + "/" + query + ".log",
                    os.system("(stretcher -asequence " + hit_output + " -sprotein1 -bsequence " +
                              query_output + " -auto -aformat clustal -stdout > " + clustal_output + ") 2> " +log_output)

                    png_output = "results/discarded/alignments/" + subject + "/" + query + ".png"
                    os.system("(jalview -nodisplay -props " + params[0] + " -colour clustal -open " +
                              clustal_output + " -png " + png_output + ") 2> " + log_output)
        except Exception as ex:
            with open("log/results/discarded/alignments.log", "w") as log_writer:
                log_writer.write(str(ex))

        os.system("touch finished.txt")
