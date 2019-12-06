import os
import traceback
from Bio import SeqIO
import multiprocessing as mp
from functools import partial


def overlap(start1, end1, lenght1, start2, end2, lenght2):
    if((start1 <= end1 and start2 >= end2) or
       (start1 >= end1 and start2 <= end2)):
        return 0

    if(start1 <= start2 and end1 >= end2):
        return 100

    if(start1 >= start2 and end1 <= end2):
        return 100

    if(start1 > end1):
        return ((min(start1, start2) - max(end1, end2)) / min(lenght1, lenght2)) * 100
    else:
        return ((min(end1, end2) - max(start1, start2)) / min(lenght1, lenght2)) * 100


def joinResultsMultiprocessing(sc, id_threshold, sim_threshold, ol_percentage):
    retained_list = []
    discarded_list = []
    result_dic = {}
    index = -1
    index_remove_list = []
    if(os.stat(sc).st_size != 0):
        subject = sc.split("/")[-2]
        subjects = SeqIO.index("data/subjects/" + subject + ".fna", "fasta")
        with open(sc, "r") as result_reader:
            content = result_reader.readlines()
            query = None
            for line in content:
                stripped_line = line.strip()
                if(stripped_line):
                    if(stripped_line.startswith("#")):
                        for contig in result_dic:
                            for index,det in result_dic[contig].items():
                                if(float(det[5]) >= id_threshold or float(det[6]) >= sim_threshold):
                                    retained_list.append("\t".join(det) + "\t" + str(index))
                                else:
                                    discarded_list.append("\t".join(det) + "\t" + str(index))

                        query = stripped_line[1:].strip()
                        result_dic.clear()
                    else:
                        header = stripped_line.split("\t")[0].strip()
                        contig = header.split("_cstart:")[0].strip()
                        contig_start = int(header.split("_cstart:")[-1].split("_cend:")[0].strip())
                        hit_start = int(header.split("::hstart=")[-1].split("::hend=")[0].strip()) + contig_start
                        hit_end = int(header.split("::hend=")[-1].split("::query=")[0].strip()) + contig_start
                        identity = float(stripped_line.split("\t")[1].strip())
                        similarity = float(stripped_line.split("\t")[2].strip())
                        hit_seq = stripped_line.split("\t")[3].strip()
                        query_seq = stripped_line.split("\t")[4].strip()
                        if(not contig in result_dic):
                            index += 1
                            result_dic[contig] = {}
                            result_dic[contig][index] = [subject, query, contig, str(hit_start), str(hit_end), str(identity), str(similarity), hit_seq, query_seq]
                        else:
                            no_overlap = True
                            for index_key,det in result_dic[contig].items():
                                if(overlap(hit_start, hit_end, len(hit_seq), int(det[3]), int(det[4]), len(det[7])) >= ol_percentage):
                                    no_overlap = False
                                    if((hit_seq.startswith("M") and det[7].startswith("M")) or
                                       not (hit_seq.startswith("M") or det[7].startswith("M"))):
                                        if(identity > float(det[5])):
                                            index_remove_list.append(index_key)
                                        elif(identity == float(det[5]) and similarity > float(det[6])):
                                            index_remove_list.append(index_key)
                                    elif(hit_seq.startswith("M")):
                                        index_remove_list.append(index_key)

                            if(no_overlap):
                                index += 1
                                result_dic[contig][index] = [subject, query, contig, str(hit_start), str(hit_end), str(identity), str(similarity), hit_seq, query_seq]
                            elif(len(index_remove_list)):
                                for index_key in index_remove_list:
                                    result_dic[contig].pop(index_key, None)

                                index += 1
                                result_dic[contig][index] = [subject, query, contig, str(hit_start), str(hit_end), str(identity), str(similarity), hit_seq, query_seq]

                            del index_remove_list[:]

    result_dic.clear()
    return[retained_list, discarded_list]


rule Remove_Duplicates:
    input:
        expand("scores/{subject}/{query}.sim", subject=config["subjects"], query=config["queries"])
    output:
        "results/retained/proda.tsv",
        "results/discarded/proda.tsv"
    params:
        config["id_threshold"],
        config["sim_threshold"],
        config["overlap_percentage"]
    threads: config["threads"]
    log:
        "log/results/results.log"
    run:
        try:
            pool = mp.Pool(processes=threads)
            pool_map = partial(joinResultsMultiprocessing, id_threshold=params[0], sim_threshold=params[1], ol_percentage=params[2])
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
                                          "\tHit sequence\tQuery sequence\tIndex\n" + "\n".join(retained_joined_results))

                with open(output[1], "w") as discarded_writer:
                    discarded_writer.write("Subject\tQuery\tContig\tStart\tEnd\tIdentity\tSimilarity"
                                           "\tHit sequence\tQuery sequence\tIndex\n" + "\n".join(discarded_joined_results))

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
