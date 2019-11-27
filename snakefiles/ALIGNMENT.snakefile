import os
from Bio import SeqIO
import traceback
import multiprocessing as mp
from functools import partial

manager = mp.Manager()
current_progress = manager.Value("i", 0.0)

def stretcherAlignmentMultiprocessing(query, targets, output, log):
    temp_query = output[0].replace(".sim", "_" + query.id + ".faa")
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + query.id + "\n" + str(query.seq))

    scores = []
    query_found = False
    for target in targets:
        query_id = target.id.split("::query=")[-1]
        if(query.id == query_id):
            if(not query_found):
                scores.append("#" + query.id)
                query_found = True

            temp_target = output[0].replace(".sim", "_target.faa")
            with open(temp_target, "w") as target_writer:
                target_writer.write(">" + target.id + "\n" + str(target.seq))

            temp_output = output[0].replace(".sim", "_" + query.id + "_target.st")
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
    current_progress.value += 1
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
                queries = list(SeqIO.parse(open(input[0]), "fasta"))
                targets = list(SeqIO.parse(open(input[1]), "fasta"))
                al_threads = threads
                if(len(queries) < al_threads):
                    al_threads = len(queries)

                pool = mp.Pool(processes=al_threads)
                pool_map = partial(stretcherAlignmentMultiprocessing, targets=targets, output=output, log=log[0])
                al_mp_results = pool.map_async(pool_map, queries)
                pool.close()
                while(current_progress.value != len(queries)):
                    if(al_mp_results.ready()):
                        pool.terminate()
                        raise Exception(al_mp_results.get()[0])

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
