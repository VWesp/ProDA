import os
from Bio import SeqIO
import traceback
import multiprocessing as mp
from functools import partial


def stretcherAlignmentMultiprocessing(query, targets, log):
    temp_query = output[0].replace(".sim", "_" + query.id + ".faa")
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + query.id + "\n" + str(query.seq))

    scores = []
    query_found = False
    for target in targets:
        query_id = target.id.split("::query=")[-1]
        if(query.id == query_id):
            if(not query_found):
                scores = ["#" + query.id]
                query_found = True

            temp_target = output[0].replace(".sim", "_" + target.id + ".faa")
            with open(temp_target, "w") as target_writer:
                target_writer.write(">" + target.id + "\n" + str(target.seq))

            temp_output = output[0].replace(".sim", "_" + query.id + "_" + target.id + ".st")
            os.system("(stretcher -asequence " + temp_query + " -sprotein1 -bsequence " +
                      temp_target + " -auto -stdout > " + temp_output + ") 2> " + log)
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
    log:
        "log/scores/{subject}/{query}.log"
    run:
        try:
            if __name__ == '__main__':
                if(os.stat(input[1]).st_size != 0):
                    queries = SeqIO.parse(open(input[0]), "fasta")
                    targets = SeqIO.parse(open(input[1]), "fasta")
                    al_threads = threads[0]
                    if(len(queries) < al_threads):
                        al_threads = len(queries)

                    pool = mp.Pool(processes=al_threads)
                    pool_map = partial(stretcherAlignmentMultiprocessing, targets=targets, log=log[0])
                    al_mp_results = pool.map_async(pool_map, queries)
                    pool.close()
                    manager = mp.Manager()
                    current_progress = manager.Value("i", 0.0)
                    while(current_progress.value != len(queries)):
                        if(al_mp_results.ready()):
                            pool.terminate()
                            raise Exception()

                    pool.join()
                    al_joined_results = [r for result in al_mp_results.get() for r in result]
                    with open(output[0], "w") as score_writer:
                        score_writer.write("\n\n".join(al_joined_results))

                    del scores[:]
                    del al_mp_results.get()[:]

                if(not os.path.exists(output[0])):
                    os.system("touch " + output[0])
        except:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()))
