import os
from Bio import SeqIO

rule exonerate:
    input:
        queries = expand("data/queries/{query}.faa", query=config["queries"]),
        blast_results = expand("matches/{subject}/{query}.fna", subject=config["subjects"], query=config["queries"])
    output:
        "exonerate/{subject}/{query}.faa"
    params:
        per=config["exonerate_percentage"]
    log:
        "log/exonerate/{subject}/{query}.log"
    run:
        for result in input.blast_results:
            result_subject_name = result.split("/")[-2]
            result_query_name = result.split("/")[-1].split(".fna")[0]
            for query in input.queries:
                query_name = query.split("/")[-1].split(".faa")[0]
                if(result_query_name == query_name):
                    exonerate_output = "exonerate/" + result_subject_name + "/" + query_name + ".ryo"
                    exonerate_log = "log/exonerate/" + result_subject_name + "/" + query_name + ".log"
                    if(not os.path.exists("exonerate/" + result_subject_name)):
                        os.makedirs("exonerate/" + result_subject_name)

                    if(not os.path.exists("log/exonerate/" + result_subject_name)):
                        os.makedirs("log/exonerate/" + result_subject_name)

                    os.system("(exonerate --model protein2genome --targettype dna --querytype protein "
                              "--ryo '>%ti::orientation=%g::score=%s::similarity=%ps::query=%qi\n%tcs' "
                              "--showalignment no --showvulgar no --refine region --percent " + str(params[0]) + " "
                              "--query " + query + " --target " + result +
                              " > " + exonerate_output + ") 2> " + exonerate_log)
                    exonerate_fasta_output = "exonerate/" + result_subject_name + "/" + query_name + ".fna"
                    os.system("(tail -n +4 " + exonerate_output + " | head -n -1) > " + exonerate_fasta_output)
                    translated_fastas = []
                    fastas = SeqIO.parse(open(exonerate_fasta_output), "fasta")
                    for fasta in fastas:
                        translated_fastas.append(">" + fasta.id + "\n" + str(fasta.seq.translate()))

                    translated_exonerate_fasta_output = "exonerate/" + result_subject_name + "/" + query_name + ".faa"
                    with open(translated_exonerate_fasta_output, "w") as translated_writer:
                        translated_writer.write("\n".join(translated_fastas))

                    del translated_fastas[:]
