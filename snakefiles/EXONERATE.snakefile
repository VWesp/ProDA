import os
from Bio.Seq import Seq
import traceback
import multiprocessing as mp
from functools import partial


manager = mp.Manager()
index = manager.Value("i", 0)
lock = manager.Lock()

def exonerateSearchMultiprocessing(match, input, blosum, percent, output, log):
    queries = SeqIO.index(input, "fasta")
    local_index = None
    with lock:
        index.value += 1
        local_index = index.value

    temp_target = output[0].replace(".gff", "_" + str(local_index) + "_target.fna")
    with open(temp_target, "w") as target_writer:
        target_writer.write(">" + match.id + "\n" + str(match.seq))

    query = match.id.split("_query:")[-1]
    temp_query = output[0].replace(".gff", "_" + str(local_index) + "_query.faa")
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + queries[query].id + "\n" + str(queries[query].seq))

    temp_output_gff = output[0].replace(".gff", "_" + str(local_index) + "_output.gff")
    os.system("(exonerate --model protein2genome --targettype dna --querytype protein "
              "--showtargetgff --showalignment no --showvulgar no "
              "--refine region --proteinsubmat blosum/" + str(blosum) + ".txt --percent " + str(percent) +
              " --query " + temp_query + " --target " + temp_target +
              " > " + temp_output_gff + ") 2> " + log)
    temp_output_pgff = output[0].replace(".gff", "_" + str(local_index) + "_output.pgff")
    os.system("(tail -n +4 " + temp_output_gff + " | head -n -1) > " + temp_output_pgff)
    ryo_results = None
    with open(temp_output_pgff, "r") as output_reader:
        ryo_results = output_reader.readlines()

    os.remove(temp_target)
    os.remove(temp_output_gff)
    os.remove(temp_output_pgff)
    os.remove(temp_query)
    return ryo_results


rule Build_Exonerate_Alignment:
    input:
        "data/queries/{query}.faa",
        "matches/{subject}/{query}.fna"
    output:
        "alignment/exonerate/{subject}/{query}.gff"
    params:
        config["blosum"],
        config["exonerate_percentage"]
    threads: config["threads"]
    log:
        "log/alignment/exonerate/{subject}/{query}.log"
    run:
        try:
            if(os.stat(input[1]).st_size != 0):
                subjects = SeqIO.index(input[0], "fasta")
                matches = list(SeqIO.parse(input[1], "fasta"))
                pool = mp.Pool(processes=threads)
                pool_map = partial(exonerateSearchMultiprocessing, input=input[0], blosum=params[0],
                                   percent=params[1], output=output, log=log[0])
                ryo_mp_results = pool.map_async(pool_map, matches)
                pool.close()
                pool.join()
                ryo_joined_results = [r.strip() for result in ryo_mp_results.get() for r in result]
                with open(output[0], "w") as ryo_writer:
                    ryo_writer.write("\n".join(ryo_joined_results))

                del ryo_joined_results[:]
                del ryo_mp_results.get()[:]

            if(not os.path.exists(output[0])):
                os.system("touch " + output[0])
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()) + "\n" + str(ex))
