import os
from Bio import SeqIO
import traceback
import multiprocessing as mp
from functools import partial

def stretcherAlignmentMultiprocessing(query, targets, output, log):
    temp_query = output[0].replace(".sim", "_" + query.id + ".faa")
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + query.id + "\n" + str(query.seq))

    scores = []
    query_found = False
    for target in SeqIO.parse(open(targets), "fasta"):
        query_id = target.id.split("::query=")[-1]
        if(query.id == query_id):
            if(not query_found):
                scores.append("#" + query.id)
                query_found = True

            temp_target = output[0].replace(".sim", "_" + query.id + "_target.fna")
            with open(temp_target, "w") as target_writer:
                target_writer.write(">" + target.id + "\n" + str(target.seq))

            temp_output = output[0].replace(".sim", "_" + query.id + "_output.st")
            os.system("(stretcher -asequence " + temp_query + " -sprotein1 -bsequence " +
                      temp_target + " -auto -stdout > " + temp_output + ") 2> " + log)
            with open(temp_output, "r") as output_reader:
                content = output_reader.readlines()
                idenity = None
                for line in content:
                    if(line.startswith("# Identity")):
                        idenity = line.split("(")[-1][:-3]
                    if(line.startswith("# Similarity")):
                        similarity = line.split("(")[-1][:-3]
                        scores.append(target.id +"\t" + idenity + "\t" + similarity + "\t" + str(target.seq) + "\t" + str(query.seq))
                        break

                scores.append(output_reader.read())

            os.remove(temp_target)
            os.remove(temp_output)

    os.remove(temp_query)
    return scores


rule Retrieve_Sequence_Similarities:
    input:
        "data/queries/{query}.faa",
        "cd_hit/{subject}/{query}_merged.faa"
    output:
        "scores/{subject}/{query}.sim"
    threads: config["threads"]
    log:
        "log/scores/{subject}/{query}.log"
    run:
        try:
            if(os.stat(input[1]).st_size != 0):
                subject = input[1].split("/")[-2]
                query = input[1].split("/")[-1].split(".fna")[0]
                queries = list(SeqIO.parse(open(input[0]), "fasta"))
                pool = mp.Pool(processes=threads)
                pool_map = partial(stretcherAlignmentMultiprocessing, targets=input[1], output=output, log=log[0])
                al_mp_results = pool.map_async(pool_map, queries)
                pool.close()
                pool.join()
                al_joined_results = [r.strip() for result in al_mp_results.get() for r in result]
                with open(output[0], "w") as score_writer:
                    score_writer.write("\n\n".join(al_joined_results))

                del al_joined_results[:]
                del al_mp_results.get()[:]

            if(not os.path.exists(output[0])):
                os.system("touch " + output[0])
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()) + "\n" + str(ex))
