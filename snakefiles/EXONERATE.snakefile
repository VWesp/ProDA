import os
from Bio.Seq import Seq
import traceback
import multiprocessing as mp
from functools import partial


manager = mp.Manager()
index = manager.Value("i", 0)
lock = manager.Lock()

def exonerateSearchMultiprocessing(query, matches, blosum, percent, output, log):
    local_index = None
    with lock:
        index.value += 1
        local_index = index.value

    temp_query = output[0].replace(".faa", "_" + str(local_index) + "_query.faa")
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + query.id + "\n" + str(query.seq))

    ryo_results = []
    for match in SeqIO.parse(matches, "fasta"):
        query_id = match.id.split("_query:")[-1]

        if(query.id == query_id):
            with lock:
                index.value += 1
                local_index = index.value

            temp_target = output[0].replace(".faa", "_" + str(local_index) + "_target.fna")
            with open(temp_target, "w") as target_writer:
                target_writer.write(">" + match.id + "\n" + str(match.seq))

            temp_output_ryo = output[0].replace(".faa", "_" + str(local_index) + "_output.ryo")
            os.system("(exonerate --model protein2genome --targettype dna --querytype protein "
                      "--ryo '>%ti::hstart=%tcb::hend=%tce::query=%qi\n%tcs' --showalignment no --showvulgar no "
                      "--refine region --proteinsubmat blosum/" + str(blosum) + ".txt --percent " + str(percent) +
                      " --query " + temp_query + " --target " + temp_target +
                      " > " + temp_output_ryo + ") 2> " + log)
            temp_output_fna = output[0].replace(".faa", "_" + str(local_index) + "_output.fna")
            os.system("(tail -n +4 " + temp_output_ryo + " | head -n -1) > " + temp_output_fna)
            with open(temp_output_fna, "r") as output_reader:
                ryo_results.append(output_reader.read())

            os.remove(temp_target)
            os.remove(temp_output_ryo)
            os.remove(temp_output_fna)

    os.remove(temp_query)
    return list(filter(None, ryo_results))


rule Build_Exonerate_Alignment:
    input:
        "data/queries/{query}.faa",
        "matches/{subject}/{query}.fna"
    output:
        "exonerate/{subject}/{query}.faa"
    params:
        config["blosum"],
        config["exonerate_percentage"]
    threads: config["threads"]
    log:
        "log/exonerate/{subject}/{query}.log"
    run:
        try:
            if(os.stat(input[1]).st_size != 0):
                queries = list(SeqIO.parse(input[0], "fasta"))
                pool = mp.Pool(processes=threads)
                pool_map = partial(exonerateSearchMultiprocessing, matches=input[1], blosum=params[0],
                                   percent=params[1], output=output, log=log[0])
                ryo_mp_results = pool.map_async(pool_map, queries)
                pool.close()
                pool.join()
                ryo_joined_results = [r.strip() for result in ryo_mp_results.get() for r in result]
                translated_fastas = []
                for result in ryo_joined_results:
                    hits = result.split(">")
                    for hit in hits:
                        lines = hit.split("\n")
                        aa_sequence = str(Seq("".join(lines[1:])).translate())
                        if(not "*" in aa_sequence):
                            translated_fastas.append(">" + lines[0] + "\n" + aa_sequence)

                with open(output[0], "w") as ryo_writer:
                    ryo_writer.write("\n".join(translated_fastas))

                del ryo_joined_results[:]
                del ryo_mp_results.get()[:]
                del translated_fastas[:]

            if(not os.path.exists(output[0])):
                os.system("touch " + output[0])
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()) + "\n" + str(ex))
