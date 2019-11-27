import os
from Bio import SeqIO
import traceback
import multiprocessing as mp
from functools import partial


def spalnSearchMultiprocessing(query, matches, pam, log):
    temp_query = output[1].replace(".faa", "_" + query.id + ".faa")
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + query.id + "\n" + str(query.seq))

    sp_results = []
    for match in matches:
        query_id = match.id.split("_query:")[-1]
        if(query.id == query_id):
            temp_target = output[1].replace(".faa", "_" + match.id + ".fna")
            with open(temp_target, "w") as target_target:
                target_target.write(">" + match.id + "\n" + str(match.seq))

            temp_output = output[0].replace(".sp", "_" + query.id + "_" + match.id + ".sp")
            os.system("(spaln -M -Q3 -O6 -S3 -yp" + str(pam) + " -yq" + str(pam) +
                      " -o" + temp_output + " " + temp_target + " " + temp_query + ") 2> " + log)
            with open(temp_output, "r") as ouput_reader:
                sp_results.append(ouput_reader.read())

            os.remove(temp_target)
            os.remove(temp_output)

    os.remove(temp_query)
    current_progress.value += 1
    return sp_results


rule spaln:
    input:
        "data/queries/{query}.faa",
        "matches/{subject}/{query}.fna"
    output:
        "spaln/{subject}/{query}.sp",
        "spaln/{subject}/{query}.faa"
    params:
        config["pam"]
    threads: config["threads"]
    log:
        "log/spaln/{subject}/{query}.log"
    run:
        try:
            if __name__ == '__main__':
                if(os.stat(input[1]).st_size != 0):
                    queries = list(SeqIO.parse(open(input[0]), "fasta"))
                    matches = SeqIO.parse(open(input[1]), "fasta")
                    sp_threads = threads[0]
                    if(len(queries) < sp_threads):
                        sp_threads = len(queries)

                    pool = mp.Pool(processes=sp_threads)
                    pool_map = partial(spalnSearchMultiprocessing, matches=matches, pam=params[0],
                                       log=log[0])
                    sp_mp_results = pool.map_async(pool_map, queries)
                    pool.close()
                    manager = mp.Manager()
                    current_progress = manager.Value("i", 0.0)
                    while(current_progress.value != len(queries)):
                        if(sp_mp_results.ready()):
                            pool.terminate()
                            raise Exception()

                    pool.join()
                    sp_joined_results = [r for result in sp_mp_results.get() for r in result]
                    with open(output[0], "w") as spaln_writer:
                        spaln_writer.write("\n".join(sp_joined_results))

                    del sp_joined_results[:]
                    del sp_mp_results.get()[:]
                    translated_fastas = []
                    with open(output[0], "r") as spaln_reader:
                        content = spaln_reader.readlines()
                        sequence = []
                        current_header = None
                        current_query = None
                        for line in content:
                            if(line.startswith(">")):
                                if(current_header != None):
                                    translated_fastas.append(">" + current_header + "\n" + str(Seq("".join(sequence)).translate()))

                                current_header = line.split(" ")[1] + "::query=" + line.split(" ")[10]
                                del sequence[:]
                            elif(line.strip().isalpha()):
                                sequence.append(line.strip())

                        translated_fastas.append(">" + current_header + "\n" + str(Seq("".join(sequence)).translate()))
                        del sequence[:]

                    with open(output[1], "w") as translated_writer:
                        translated_writer.write("\n".join(translated_fastas))

                    del translated_fastas[:]

                if(not os.path.exists(output[0])):
                    os.system("touch " + output[0])

                if(not os.path.exists(output[1])):
                    os.system("touch " + output[1])
        except:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()))
