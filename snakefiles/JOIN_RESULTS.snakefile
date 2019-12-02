import os
import traceback
from Bio import SeqIO
import multiprocessing as mp
from functools import partial


def joinResultsMultiprocessing(tsv, id_threshold, sim_threshold):
    retained_list = []
    discarded_list = []
    if(os.stat(tsv).st_size != 0):
        subject = tsv.split("/")[-2]
        subjects = SeqIO.index("data/subjects/" + subject + ".fna", "fasta")
        with open(tsv, "r") as result_reader:
            content = result_reader.readlines()
            for line in content:
                stripped_line = line.strip()
                if(stripped_line):
                    query = stripped_line.split("\t")[0].strip()
                    header = stripped_line.split("\t")[1].strip()
                    contig = header.split("_cstart:")[0].strip()
                    contig_start = int(header.split("_cstart:")[-1].split("_cend:")[0].strip())
                    hit_start = int(header.split("::hstart=")[-1].split("::hend=")[0].strip()) + contig_start
                    hit_end = int(header.split("::hend=")[-1].strip()) + contig_start
                    idenity = stripped_line.split("\t")[2].strip()
                    similarity = stripped_line.split("\t")[3].strip()
                    hit_seq = stripped_line.split("\t")[4].strip()
                    query_seq = stripped_line.split("\t")[5].strip()
                    result = subject + "\t" + query + "\t" + contig + "\t" + str(hit_start) + "\t" + str(hit_end) + "\t" + idenity + "\t" + similarity + "\t" + hit_seq + "\t" + query_seq
                    if(float(idenity) >= id_threshold and float(similarity) >= sim_threshold):
                        retained_list.append(result)
                    else:
                        discarded_list.append(result)

    return[retained_list, discarded_list]


rule Join_ProDA_Results:
    input:
        expand("best_hit/{subject}/{query}.tsv", subject=config["subjects"], query=config["queries"])
    output:
        "results/retained/proda.tsv",
        "results/discarded/proda.tsv"
    params:
        config["id_threshold"],
        config["sim_threshold"]
    threads: config["threads"]
    log:
        "log/results/results.log"
    run:
        try:
            pool = mp.Pool(processes=threads)
            pool_map = partial(joinResultsMultiprocessing, id_threshold=params[0], sim_threshold=params[1])
            jr_mp_results = pool.map_async(pool_map, input)
            pool.close()
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
                    retained_writer.write("Subject\tQuery\tContig\tStart\tEnd\tIdentity\tSimilarity"
                                          "\tHit sequence\tQuery sequence\n" + "\n".join(retained_joined_results))

                with open(output[1], "w") as discarded_writer:
                    discarded_writer.write("Subject\tQuery\tContig\tStart\tEnd\tIdentity\tSimilarity"
                                           "\tHit sequence\tQuery sequence\n" + "\n".join(discarded_joined_results))

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
