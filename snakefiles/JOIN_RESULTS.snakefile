import os

rule join_results:
    input:
        expand("best_hit/{subject}/{query}.tsv", subject=config["subjects"], query=config["queries"])
    output:
        "results/retained/proda.tsv",
        "results/discarded/proda.tsv"
    params:
        threshold=config["sim_threshold"]
    log:
        "log/results/results.log"
    run:
        try:
            retained_list = []
            discarded_list = []
            for tsv in input:
                if(os.stat(tsv).st_size != 0):
                    subject = tsv.split("/")[-2]
                    with open(tsv, "r") as result_reader:
                        content = result_reader.readlines()
                        for line in content:
                            query = line.split("\t")[0]
                            contig = line.split("\t")[1]
                            similarity = line.split("\t")[2]
                            hit_seq = line.split("\t")[3]
                            query_seq = line.split("\t")[4]
                            result = subject + "\t" + query + "\t" + contig + "\t" + similarity + "\t" + hit_seq + "\t" + query_seq
                            if(float(similarity) >= params[0]):
                                retained_list.append(result)
                            else:
                                discarded_list.append(result)

            with open(output[0], "w") as retained_writer:
                retained_writer.write("Subject\tQuery\tContig\tSimilarity\tHit sequence\tQuery sequence\n" + "\n".join(retained_list))

            with open(output[1], "w") as discarded_writer:
                discarded_writer.write("Subject\tQuery\tContig\tSimilarity\tSequence\tQuery sequence\n" + "\n".join(discarded_list))

            del retained_list[:]
            del discarded_list[:]
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(ex))
