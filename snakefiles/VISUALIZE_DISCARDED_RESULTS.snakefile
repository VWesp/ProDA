import os
import traceback
import multiprocessing as mp


def visualizeDiscardedResultsMultiprocessing(line, properties):
    stripped_line = line.strip()
    if(stripped_line):
        subject = stripped_line.split("\t")[0]
        query = stripped_line.split("\t")[1]
        hit_sequence = stripped_line.split("\t")[7]
        query_sequence = stripped_line.split("\t")[8]
        if(not os.path.exists("results/discarded/alignments/" + subject)):
            os.makedirs("results/discarded/alignments/" + subject)

        if(not os.path.exists("log/results/discarded/alignments/" + subject)):
            os.makedirs("log/results/discarded/alignments/" + subject)

        hit_output =  "results/discarded/alignments/" + subject + "/" + query + "_hit.faa"
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
        os.system("(jalview -nodisplay -props " + properties + " -colour clustal -open " +
                  clustal_output + " -png " + png_output + ") 2> " + log_output)

        os.remove(hit_output)
        os.remove(query_output)


rule Visualize_Discarded_Results:
    input:
        "results/discarded/proda.tsv"
    output:
        temp("discarded_finished.txt")
    params:
        props=config["properties"]
    threads: config["threads"]
    run:
        try:
            if(os.stat(input[0]).st_size != 0):
                with open(input[0], "r") as discarded_reader:
                    lines = discarded_reader.read().split("\n")[1:]
                    pool = mp.Pool(processes=threads)
                    pool_map = partial(visualizeDiscardedResultsMultiprocessing, properties=params[0])
                    pool.map_async(pool_map, lines)
                    pool.close()
                    pool.join()
        except:
            print("\033[1;31;mError: " + str(traceback.format_exc()) + "\nSee log file: log/results/discarded/alignments.log")
            if(not os.path.exists("log/results/discarded")):
                os.makedirs("log/results/discarded")

            with open("log/results/discarded/alignments.log", "w") as log_writer:
                log_writer.write(str(traceback.format_exc()))

        os.system("touch discarded_finished.txt")
