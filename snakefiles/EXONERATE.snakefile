import os
from Bio import SeqIO
import traceback
import multiprocessing as mp
from functools import partial


def exonerateSearchMultiprocessing(query, matches, blosum, percent, log):
    temp_query = output[2].replace(".faa", "_" + query.id + ".faa")
    with open(temp_query, "w") as query_writer:
        query_writer.write(">" + query.id + "\n" + str(query.seq))

    ryo_results = []
    for match in matches:
        query_id = match.id.split("_query:")[-1]
        if(query.id == query_id):
            temp_target = output[2].replace(".faa", "_" + match.id + ".faa")
            with open(temp_target, "w") as target_target:
                target_target.write(">" + match.id + "\n" + str(match.seq))

            temp_output = output[0].replace(".ryo", "_" + query.id + "_" + match.id + ".ryo")
            os.system("(exonerate --model protein2genome --targettype dna --querytype protein "
                      "--ryo '>%ti::query=%qi\n%tcs' --showalignment no --showvulgar no "
                      "--refine region --proteinsubmat blosum/" + str(blosum) + ".txt --percent " + str(percent) +
                      " --query " + temp_query + " --target " + temp_target +
                      " > " + temp_output + ") 2> " + log)

            with open(temp_output, "r") as ouput_reader:
                ryo_results.append(ouput_reader.read())

            os.remove(temp_target)
            os.remove(temp_output)

    os.remove(temp_query)
    current_progress.value += 1
    return ryo_results


rule Build_Exonerate_Alignment:
    input:
        "data/queries/{query}.faa",
        "matches/{subject}/{query}.fna"
    output:
        "exonerate/{subject}/{query}.ryo",
        "exonerate/{subject}/{query}.faa",
        temp("exonerate/{subject}/{query}.fna")
    params:
        bls=config["blosum"],
        per=config["exonerate_percentage"]
    threads: config["threads"]
    log:
        "log/exonerate/{subject}/{query}.log"
    run:
        try:
            if __name__ == '__main__':
                if(os.stat(input[1]).st_size != 0):
                    queries = list(SeqIO.parse(open(input[0]), "fasta"))
                    matches = SeqIO.parse(open(input[1]), "fasta")
                    ex_threads = threads[0]
                    if(len(queries) < ex_threads):
                        ex_threads = len(queries)

                    pool = mp.Pool(processes=ex_threads)
                    pool_map = partial(exonerateSearchMultiprocessing, matches=matches, blosum=params[0],
                                       percent=params[1], log=log[0])
                    ryo_mp_results = pool.map_async(pool_map, queries)
                    pool.close()
                    manager = mp.Manager()
                    current_progress = manager.Value("i", 0.0)
                    while(current_progress.value != len(queries)):
                        if(ryo_mp_results.ready()):
                            pool.terminate()
                            raise Exception()

                    pool.join()
                    ryo_joined_results = [r for result in ryo_mp_results.get() for r in result]
                    with open(output[0], "w") as ryo_writer:
                        ryo_writer.write("\n".join(ryo_joined_results))

                    del ryo_joined_results[:]

                    #ryo: ::orientation=%g::score=%s::similarity=%ps
                    '''os.system("(exonerate --model protein2genome --targettype dna --querytype protein "
                              "--ryo '>%ti::query=%qi\n%tcs' --showalignment no --showvulgar no "
                              "--refine region --proteinsubmat blosum/" + str(params[0]) + ".txt --percent " + str(params[1]) + " "
                              "--query " + input[0] + " --target " + input[1] +
                              " > " + output[0] + ") 2> " + log[0])'''

                    os.system("(tail -n +4 " + output[0] + " | head -n -1) > " + output[2])
                    fastas = SeqIO.parse(open(output[2]), "fasta")
                    translated_fastas = []
                    for fasta in fastas:
                        translated_fastas.append(">" + fasta.id + "\n" + str(fasta.seq.translate()))

                    with open(output[1], "w") as translated_writer:
                        translated_writer.write("\n".join(translated_fastas))

                    del translated_fastas[:]

                if(not os.path.exists(output[0])):
                    os.system("touch " + output[0])

                if(not os.path.exists(output[1])):
                    os.system("touch " + output[1])

                if(not os.path.exists(output[2])):
                    os.system("touch " + output[2])
        except:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()))
