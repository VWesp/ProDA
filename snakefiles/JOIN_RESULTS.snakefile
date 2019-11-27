import os
import traceback
import multiprocessing as mp


def joinResultsMultiprocessing(tsv):
    retained = []
    discarded = []
    if(os.stat(tsv).st_size != 0):
        subject = tsv.split("/")[-2]
        with open(tsv, "r") as result_reader:
            content = result_reader.readlines()
            for line in content:
                query = line.split("\t")[0].strip()
                contig = line.split("\t")[1].strip()
                similarity = line.split("\t")[2].strip()
                hit_seq = line.split("\t")[3].strip()
                query_seq = line.split("\t")[4].strip()
                result = subject + "\t" + query + "\t" + contig + "\t" + similarity + "\t" + hit_seq + "\t" + query_seq
                if(float(similarity) >= params[0]):
                    retained_list.append(result)
                else:
                    discarded_list.append(result)

    current_progress.value += 1
    return[retained, discarded]


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
            if __name__ == '__main__':
                jr_threads = threads[0]
                if(len(input) < jr_threads):
                    jr_threads = len(input)

                pool = mp.Pool(processes=jr_threads)
                jr_mp_results = pool.map_async(stretcherAlignmentMultiprocessing, input)
                pool.close()
                manager = mp.Manager()
                current_progress = manager.Value("i", 0.0)
                while(current_progress.value != len(queries)):
                    if(jr_mp_results.ready()):
                        pool.terminate()
                        raise Exception()

                pool.join()
                retained_joined_results = [r for result in jr_mp_results.get() for r[0] in result]
                discarded_joined_results = [r for result in jr_mp_results.get() for r[1] in result]
                if(len(retained_joined_results) or len(discarded_joined_results)):
                    with open(output[0], "w") as retained_writer:
                        retained_writer.write("Subject\tQuery\tContig\tSimilarity\tHit sequence\tQuery sequence\n" + "\n".join(retained_joined_results))

                    with open(output[1], "w") as discarded_writer:
                        discarded_writer.write("Subject\tQuery\tContig\tSimilarity\tSequence\tQuery sequence\n" + "\n".join(discarded_joined_results))

                    del retained_joined_results[:]
                    del discarded_joined_results[:]
                    del jr_mp_results.get()[:]

                if(not os.path.exists(output[0])):
                    os.system("touch " + output[0])

                if(not os.path.exists(output[1])):
                    os.system("touch " + output[1])
        except:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()))
