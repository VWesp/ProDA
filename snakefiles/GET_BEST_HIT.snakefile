import os

rule get_best_hit:
    input:
        "scores/{subject}/{query}.sim"
    output:
        "best_hit/{subject}/{query}.tsv"
    log:
        "log/best_hit/{subject}/{query}.log"
    run:
        try:
            if(os.stat(input[0]).st_size != 0):
                similarities = {}
                with open(input[0], "r") as score_reader:
                    content = score_reader.readlines()
                    current_query = None
                    for line in content:
                        if(line.strip()):
                            if(line.startswith("#")):
                                current_query = line[1:].strip()
                                similarities[current_query] = ["", -1, "", ""]
                            else:
                                contig = line.split("\t")[0].strip()
                                similarity = float(line.split("\t")[1].strip())
                                hit_sequence = line.split("\t")[2].strip()
                                query_sequence = line.split("\t")[3].strip()
                                if(similarity > similarities[current_query][1]):
                                    similarities[current_query] = [contig, similarity, hit_sequence, query_sequence]
                best_hits = []
                if(similarities):
                    for query, hit in similarities.items():
                        contig = hit[0].split("::query=")[0]
                        best_hits.append(query + "\t" + contig + "\t" + str(hit[1]) + "\t" + hit[2] + "\t" + hit[3])

                with open(output[0], "w") as statistic_writer:
                    statistic_writer.write("\n".join(best_hits))

                similarities.clear()
                del best_hits[:]

            if(not os.path.exists(output[0])):
                os.system("touch " + output[0])
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(ex))
