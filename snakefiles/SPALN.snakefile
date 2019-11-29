import os
from Bio import SeqIO
import traceback
import multiprocessing as mp
from functools import partial

def spalnSearchMultiprocessing(query, matches, pam, output, log):
    temp_query = output[0].replace(".faa", "_" + query.id + ".faa")
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + query.id + "\n" + str(query.seq))

    sp_results = []
    for match in SeqIO.parse(open(matches), "fasta"):
        query_id = match.id.split("_query:")[-1]
        if(query.id == query_id):
            temp_target = output[0].replace(".faa", "_" + query.id + "_target.fna")
            with open(temp_target, "w") as target_writer:
                target_writer.write(">" + match.id + "\n" + str(match.seq))

            temp_output = output[0].replace(".faa", "_" + query.id + "_output.sp")
            os.system("(spaln -M -Q3 -O6 -S3 -yp" + str(pam) + " -yq" + str(pam) +
                      " -o" + temp_output + " " + temp_target + " " + temp_query + ") 2> " + log)
            with open(temp_output, "r") as output_reader:
                sp_results.append(output_reader.read())

            os.remove(temp_target)
            os.remove(temp_output)

    os.remove(temp_query)
    return sp_results


rule Build_Spaln_Alignment:
    input:
        "data/queries/{query}.faa",
        "matches/{subject}/{query}.fna"
    output:
        "spaln/{subject}/{query}.faa",
        temp("spaln/{subject}/{query}.sp")
    params:
        config["pam"]
    threads: config["threads"]
    log:
        "log/spaln/{subject}/{query}.log"
    run:
        try:
            if(os.stat(input[1]).st_size != 0):
                subject = input[1].split("/")[-2]
                query = input[1].split("/")[-1].split(".fna")[0]
                queries = list(SeqIO.parse(open(input[0]), "fasta"))
                pool = mp.Pool(processes=threads)
                pool_map = partial(spalnSearchMultiprocessing, matches=input[1], pam=params[0],
                                   output=output, log=log[0])
                sp_mp_results = pool.map_async(pool_map, queries)
                pool.close()
                pool.join()
                sp_joined_results = [r.strip() for result in sp_mp_results.get() for r in result]
                with open(output[1], "w") as spaln_writer:
                    spaln_writer.write("\n".join(sp_joined_results))

                del sp_joined_results[:]
                del sp_mp_results.get()[:]
                translated_fastas = []
                with open(output[1], "r") as spaln_reader:
                    content = spaln_reader.readlines()
                    sequence = []
                    current_header = None
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

                with open(output[0], "w") as translated_writer:
                    translated_writer.write("\n".join(translated_fastas))

                del translated_fastas[:]

            if(not os.path.exists(output[0])):
                os.system("touch " + output[0])

            if(not os.path.exists(output[1])):
                os.system("touch " + output[1])
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()) + "\n" + str(ex))
