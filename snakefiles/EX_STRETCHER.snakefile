import os
from Bio import SeqIO
import traceback
import multiprocessing as mp
from functools import partial


manager = mp.Manager()
index = manager.Value("i", 0)
lock = manager.Lock()

def exStretcherAlignmentMultiprocessing(target, input, output, log):
    queries = SeqIO.index(input, "fasta")
    local_index = None
    with lock:
        index.value += 1
        local_index = index.value

    gff_result = list(filter(None, target.split("\n")))
    if(len(gff_result)):
        temp_target = output[0].replace(".sc_gff", "_" + str(local_index) + "_target.fna")
        with open(temp_target, "w") as target_writer:
            target_writer.write(">" + gff_result[0].split("\t")[0] + "\n" + gff_result[-1].split(">pep:")[-1])

        query = gff_result[0].split("\t")[1]
        temp_query = output[0].replace(".sc_gff", "_" + str(local_index) + "_query.faa")
        with open(temp_query, "w") as query_writer:
            query_writer.write(">" + queries[query].id + "\n" + str(queries[query].seq))

        temp_output = output[0].replace(".sc_gff", "_" + str(local_index) + "_output.sc_gff")
        os.system("(stretcher -asequence " + temp_query + " -sprotein1 -bsequence " +
                  temp_target + " -auto -stdout > " + temp_output + ") 2> " + log)

        with open(temp_output, "r") as output_reader:
            content = output_reader.readlines()
            identity = None
            similarity = None
            for line in content:
                if(line.startswith("# Identity")):
                    identity = line.split("(")[-1][:-3]
                if(line.startswith("# Similarity")):
                    similarity = line.split("(")[-1][:-3]
                    break

            gff_result[0] = "#" + gff_result[0] + "\t" + identity + "\t" + similarity

        os.remove(temp_target)
        os.remove(temp_output)
        os.remove(temp_query)

    return gff_result


rule Retrieve_Exonerate_Sequence_Similarities:
    input:
        "data/queries/{query}.faa",
        "alignment/exonerate/{subject}/{query}.pr_gff"
    output:
        "scores/exonerate/{subject}/{query}.sc_gff"
    threads: config["threads"]
    log:
        "log/scores/exonerate/{subject}/{query}.log"
    run:
        try:
            if(os.stat(input[1]).st_size != 0):
                subject = input[1].split("/")[-2]
                query = input[1].split("/")[-1].split(".fna")[0]
                gff_content = None
                with open(input[1], "r") as gff_reader:
                    gff_content = gff_reader.read().split("#")

                pool = mp.Pool(processes=threads)
                pool_map = partial(exStretcherAlignmentMultiprocessing, input=input[0], output=output, log=log[0])
                al_mp_results = pool.map_async(pool_map, gff_content)
                pool.close()
                pool.join()
                al_joined_results = [r.strip() for result in al_mp_results.get() for r in result]
                with open(output[0], "w") as score_writer:
                    score_writer.write("\n".join(al_joined_results))

                del al_joined_results[:]
                del al_mp_results.get()[:]

            if(not os.path.exists(output[0])):
                os.system("touch " + output[0])
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()) + "\n" + str(ex))
