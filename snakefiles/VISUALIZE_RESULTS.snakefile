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
            if(os.stat(input[0]).st_size != 0):
                with open(input[0], "r") as proda_reader:
                    next(proda_reader)
                    content = proda_reader.readlines()
                    for line in content:
                        subject = line.split("\t")[0]
                        query = line.split("\t")[1]
                        hit_sequence = line.split("\t")[4]
                        query_sequence = line.split("\t")[5]
                        hit_output =  "results/retained/alignments/" + subject + "/" + query + "_hit.faa"
                        if(not os.path.exists("results/retained/alignments/" + subject)):
                            os.makedirs("results/retained/alignments/" + subject)

                        if(not os.path.exists("log/results/retained/alignments/" + subject)):
                            os.makedirs("log/results/retained/alignments/" + subject)

                        with open(hit_output, "w") as hit_writer:
                            hit_writer.write(">" + subject + "\n" + hit_sequence)

                        query_output =  "results/retained/alignments/" + subject + "/" + query + "_query.faa"
                        with open(query_output, "w") as query_writer:
                            query_writer.write(">" + query + "\n" + query_sequence)

                        clustal_output =  "results/retained/alignments/" + subject + "/" + query + ".clustal"
                        log_output = "log/results/retained/alignments/" + subject + "/" + query + ".log"

                        os.system("(stretcher -asequence " + hit_output + " -sprotein1 -bsequence " +
                                  query_output + " -auto -aformat clustal -stdout > " + clustal_output + ") 2> " +log_output)

                        png_output = "results/retained/alignments/" + subject + "/" + query + ".png"
                        os.system("(jalview -nodisplay -props " + params[0] + " -colour clustal -open " +
                                  clustal_output + " -png " + png_output + ") 2> " + log_output)

                        if(os.path.exists(hit_output)):
                            os.remove(hit_output)

                        if(os.path.exists(query_output)):
                            os.remove(query_output)

                        if(os.path.exists(query_output)):
                            os.remove(query_output)
        except Exception as ex:
            print("\033[1;31;mError: " + str(ex) + "\nSee log file: log/results/retained/alignments.log")
            if(not os.path.exists("log/results/retained")):
                os.makedirs("log/results/retained")

            with open("log/results/retained/alignments.log", "w") as log_writer:
                log_writer.write(str(ex))

        try:
            if(os.stat(input[1]).st_size != 0):
                with open(input[1], "r") as proda_reader:
                    next(proda_reader)
                    content = proda_reader.readlines()
                    for line in content:
                        subject = line.split("\t")[0]
                        query = line.split("\t")[1]
                        hit_sequence = line.split("\t")[4]
                        query_sequence = line.split("\t")[5]
                        hit_output =  "results/discarded/alignments/" + subject + "/" + query + "_hit.faa"
                        if(not os.path.exists("results/discarded/alignments/" + subject)):
                            os.makedirs("results/discarded/alignments/" + subject)

                        if(not os.path.exists("log/results/discarded/alignments/" + subject)):
                            os.makedirs("log/results/discarded/alignments/" + subject)

                        with open(hit_output, "w") as hit_writer:
                            hit_writer.write(">" + subject + "\n" + hit_sequence)

                        query_output =  "results/discarded/alignments/" + subject + "/" + query + "_query.faa"
                        with open(query_output, "w") as query_writer:
                            query_writer.write(">" + query + "\n" + query_sequence)

                        clustal_output =  "results/discarded/alignments/" + subject + "/" + query + ".clustal"
                        log_output = "log/results/discarded/alignments/" + subject + "/" + query + ".log"
                        os.system("(stretcher -asequence " + hit_output + " -sprotein1 -bsequence " +
                                  query_output + " -auto -aformat clustal -stdout > " + clustal_output + ") 2> " +log_output)

                        png_output = "results/discarded/alignments/" + subject + "/" + query + ".png"
                        os.system("(jalview -nodisplay -props " + params[0] + " -colour clustal -open " +
                                  clustal_output + " -png " + png_output + ") 2> " + log_output)

                        if(os.path.exists(hit_output)):
                            os.remove(hit_output)

                        if(os.path.exists(query_output)):
                            os.remove(query_output)

                        if(os.path.exists(query_output)):
                            os.remove(query_output)
        except Exception as ex:
            print("\033[1;31;mError: " + str(ex) + "\nSee log file: log/results/discarded/alignments.log")
            if(not os.path.exists("log/results/discarded")):
                os.makedirs("log/results/discarded")

            with open("log/results/discarded/alignments.log", "w") as log_writer:
                log_writer.write(str(ex))

        os.system("touch finished.txt")
