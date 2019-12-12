import os
from Bio import SeqIO
import traceback
import multiprocessing as mp
from functools import partial


manager = mp.Manager()
index = manager.Value("i", -1)
lock = manager.Lock()

def spalnSearchMultiprocessing(match, input, pam, output, log):
    queries = SeqIO.index(input, "fasta")
    local_index = None
    with lock:
        index.value += 1
        local_index = index.value

    temp_target = output[0].replace(".pr_gff", "_" + str(local_index) + "_target.fna")
    with open(temp_target, "w") as target_writer:
        target_writer.write(">" + match.id + "\n" + str(match.seq))

    query = match.id.split("_query:")[-1]
    temp_query = output[0].replace(".pr_gff", "_" + str(local_index) + "_query.faa")
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + queries[query].id + "\n" + str(queries[query].seq))

    temp_output = output[0].replace(".pr_gff", "_" + str(local_index) + "_output.gff")
    os.system("(spaln -M -Q3 -O0 -S3 -yp" + str(pam) + " -yq" + str(pam) +
              " -o" + temp_output + " " + temp_target + " " + temp_query + ") 2> " + log)

    sp_results = None
    with open(temp_output, "r") as output_reader:
        sp_results = output_reader.readlines()

    os.remove(temp_target)
    os.remove(temp_output)
    os.remove(temp_query)
    return sp_results


rule Build_Spaln_Alignment:
    input:
        "data/subjects/{subject}.fna",
        "data/queries/{query}.faa",
        "matches/{subject}/{query}.fna"
    output:
        "alignment/spaln/{subject}/{query}.pr_gff"
    params:
        config["pam"]
    threads: config["threads"]
    log:
        "log/alignment/spaln/{subject}/{query}.log"
    run:
        try:
            if(os.stat(input[2]).st_size != 0):
                subjects = SeqIO.index(input[0], "fasta")
                matches = list(SeqIO.parse(input[2], "fasta"))
                pool = mp.Pool(processes=threads)
                pool_map = partial(spalnSearchMultiprocessing, input=input[1], pam=params[0],
                                   output=output, log=log[0])
                sp_mp_results = pool.map_async(pool_map, matches)
                pool.close()
                pool.join()
                sp_joined_results = [r.strip() for result in sp_mp_results.get() for r in result]
                gff_results = readGFF(sp_joined_results, subjects)
                with open(output[0], "w") as ryo_writer:
                    ryo_writer.write("\n".join(gff_results))

                del sp_joined_results[:]
                del sp_mp_results.get()[:]
                del gff_results[:]

            if(not os.path.exists(output[0])):
                os.system("touch " + output[0])
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()) + "\n" + str(ex))
