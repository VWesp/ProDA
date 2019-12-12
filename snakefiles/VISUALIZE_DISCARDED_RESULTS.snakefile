import os
import traceback
from Bio import SeqIO
import multiprocessing as mp
from functools import partial


def visualizeDiscardedResultsMultiprocessing(hit, subject, input, properties):
    queries = SeqIO.index(input, "fasta")
    query = hit.id.split("|")[1]
    index = hit.id.split("|")[-1]
    hit_output =  "results/below_threshold/" + subject + "/alignments/" + query + "_index_" + index + "_hit.faa"
    with open(hit_output, "w") as hit_writer:
        hit_writer.write(">" + hit.id + "\n" + str(hit.seq))

    query_output =  "results/below_threshold/" + subject + "/alignments/" + query + "_index_" + index + "_query.faa"
    with open(query_output, "w") as query_writer:
        query_writer.write(">" + queries[query].id + "\n" + str(queries[query].seq))

    clustal_output =  "results/below_threshold/" + subject + "/alignments/" + query + "_index_" + index + ".clustal"
    log_output = "log/results/below_threshold/" + subject + "/alignments/" + query + "_index_" + index + ".log"
    os.system("(stretcher -asequence " + hit_output + " -sprotein1 -bsequence " +
              query_output + " -auto -aformat clustal -stdout > " + clustal_output + ") 2> " +log_output)

    png_output = "results/below_threshold/" + subject + "/alignments/" + query + "_index_" + index + ".png"
    os.system("(jalview -nodisplay -props " + properties + " -colour clustal -open " +
              clustal_output + " -png " + png_output + ") 2> " + log_output)

    os.remove(hit_output)
    os.remove(query_output)


rule Visualize_Results_Below_Threshold:
    input:
        "data/queries/{query}.faa",
        "results/below_threshold/{subject}/{query}.pep"
    output:
        temp("temp/below/{subject}/{query}.txt")
    params:
        config["properties"]
    threads: config["threads"]
    run:
        try:
            if(os.stat(input[1]).st_size != 0):
                subject = input[1].split("/")[-2]
                if(not os.path.exists("results/below_threshold/" + subject + "/alignments")):
                    os.makedirs("results/below_threshold/" + subject + "/alignments")

                if(not os.path.exists("log/results/below_threshold/" + subject + "/alignments")):
                    os.makedirs("log/results/below_threshold/" + subject + "/alignments")

                hits = list(SeqIO.parse(input[1], "fasta"))
                pool = mp.Pool(processes=threads)
                pool_map = partial(visualizeDiscardedResultsMultiprocessing, subject=subject, input=input[0], properties=params[0])
                pool.map_async(pool_map, hits)
                pool.close()
                pool.join()
        except:
            print("\033[1;31;mError: " + str(traceback.format_exc()) + "\nSee log file: log/results/below_threshold/" + subject + "/alignments/.log")
            if(not os.path.exists("log/results/below_threshold")):
                os.makedirs("log/results/below_threshold")

            with open("log/results/below_threshold/" + subject + "/alignments/.log", "w") as log_writer:
                log_writer.write(str(traceback.format_exc()))

        os.system("touch " + output[0])
