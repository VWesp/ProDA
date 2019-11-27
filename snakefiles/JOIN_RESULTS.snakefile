import os
import traceback
import multiprocessing as mp
from functools import partial

manager = mp.Manager()
current_progress = manager.Value("i", 0.0)

def joinResultsMultiprocessing(tsv, threshold):
    retained_list = []
    discarded_list = []
    if(os.stat(tsv).st_size != 0):
        subject = tsv.split("/")[-2]
        with open(tsv, "r") as result_reader:
            content = result_reader.readlines()
            for line in content:
                query = line.split("\t")[0].strip()
                contig = line.split("\t")[1].strip()
                idenity = line.split("\t")[2].strip()
                similarity = line.split("\t")[3].strip()
                hit_seq = line.split("\t")[4].strip()
                query_seq = line.split("\t")[5].strip()
                result = subject + "\t" + query + "\t" + contig + "\t" + idenity + "\t" + similarity + "\t" + hit_seq + "\t" + query_seq
                if(float(similarity) >= threshold):
                    retained_list.append(result)
                else:
                    discarded_list.append(result)

    current_progress.value += 1
    return[retained_list, discarded_list]


rule Join_ProDA_Results:
    input:
        expand("best_hit/{subject}/{query}.tsv", subject=config["subjects"], query=config["queries"])
    output:
        "results/retained/proda.tsv",
        "results/discarded/proda.tsv"
    params:
        config["sim_threshold"]
    threads: config["threads"]
    log:
        "log/results/results.log"
    run:
        try:
            jr_threads = threads
            if(len(input) < jr_threads):
                jr_threads = len(input)

            pool = mp.Pool(processes=jr_threads)
            pool_map = partial(joinResultsMultiprocessing, threshold=params[0])
            jr_mp_results = pool.map_async(pool_map, input)
            pool.close()
            while(current_progress.value != len(input)):
                if(jr_mp_results.ready()):
                    pool.terminate()
                    raise Exception(jr_mp_results.get()[0])

            pool.join()
            retained_joined_results = []
            discarded_joined_results = []
            for result in jr_mp_results.get():
                for retained in result[0]:
                    retained_joined_results.append(retained)

                for discarded in result[1]:
                    discarded_joined_results.append(discarded)

            if(len(retained_joined_results) or len(discarded_joined_results)):
                with open(output[0], "w") as retained_writer:
                    retained_writer.write("Subject\tQuery\tContig\tIdentity\tSimilarity\tHit sequence\tQuery sequence\n" + "\n".join(retained_joined_results))

                with open(output[1], "w") as discarded_writer:
                    discarded_writer.write("Subject\tQuery\tContig\tIdentity\tSimilarity\tSequence\tQuery sequence\n" + "\n".join(discarded_joined_results))

            del retained_joined_results[:]
            del discarded_joined_results[:]
            del jr_mp_results.get()[:]

            if(not os.path.exists(output[0])):
                os.system("touch " + output[0])

            if(not os.path.exists(output[1])):
                os.system("touch " + output[1])
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()) + "\n" + str(ex))
