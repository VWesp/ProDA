import os
import traceback
import multiprocessing as mp


def visualizeResultsMultiprocessing(line, properties):
    stripped_line = line.strip()
    if(stripped_line):
        subject = stripped_line.split("\t")[0]
        query = stripped_line.split("\t")[1]
        hit_sequence = stripped_line.split("\t")[5]
        query_sequence = stripped_line.split("\t")[6]
        if(not os.path.exists("results/retained/alignments/" + subject)):
            os.makedirs("results/retained/alignments/" + subject)

        if(not os.path.exists("log/results/retained/alignments/" + subject)):
            os.makedirs("log/results/retained/alignments/" + subject)

        hit_output =  "results/retained/alignments/" + subject + "/" + query + "_hit.faa"
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
        os.system("(jalview -nodisplay -props " + properties + " -colour clustal -open " +
                  clustal_output + " -png " + png_output + ") 2> " + log_output)

        os.remove(hit_output)
        os.remove(query_output)


rule Visualize_Retained_Results:
    input:
        "results/retained/proda_temp.tsv"
    output:
        temp("retained_finished.txt")
    params:
        props=config["properties"]
    threads: config["threads"]
    run:
        try:
            if(os.stat(input[0]).st_size != 0):
                with open(input[0], "r") as retained_reader:
                    lines = retained_reader.read().split("\n")[1:]
                    pool = mp.Pool(processes=threads)
                    pool_map = partial(visualizeResultsMultiprocessing, properties=params[0])
                    pool.map_async(pool_map, lines)
                    pool.close()
                    pool.join()
        except:
            print("\033[1;31;mError: " + str(traceback.format_exc()) + "\nSee log file: log/results/discarded/alignments.log")
            if(not os.path.exists("log/results/discarded")):
                os.makedirs("log/results/discarded")

            with open("log/results/discarded/alignments.log", "w") as log_writer:
                log_writer.write(str(traceback.format_exc()))

        os.system("touch retained_finished.txt")
