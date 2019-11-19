import os

rule exonerate:
    input:
        queries = expand("data/queries/{query}.faa", query=config["queries"]),
        blast_results = expand("matches/{subject}/{query}.fasta", subject=config["subjects"], query=config["queries"])
    output:
        "exonerate_results/{subject}/{query}.fasta"
    params:
        per=config["exonerate_percentage"]
    log:
        "log/exonerate_results/{subject}/{query}.log"
    run:
        for result in input.blast_results:
            result_subject_name = result.split("/")[-2]
            result_query_name = result.split("/")[-1].split(".fasta")[0]
            for query in input.queries:
                query_name = query.split("/")[-1].split(".faa")[0]
                if(result_query_name == query_name):
                    exonerate_output = "exonerate_results/" + result_subject_name + "/" + query_name + ".ryo"
                    exonerate_log = "log/exonerate_results/" + result_subject_name + "/" + query_name + ".log"
                    if(not os.path.exists("log/exonerate_results/" + result_subject_name)):
                        os.makedirs("log/exonerate_results/" + result_subject_name)

                    os.system("(exonerate --model protein2genome --targettype dna --querytype protein "
                              "--ryo '>%ti::orientation=%g::score=%s::similarity=%ps::query=%qi\n%tcs' "
                              "--showalignment no --showvulgar no --refine region --percent " + str(params[0]) + " "
                              "--query " + query + " --target " + result +
                              " > " + exonerate_output + ") 2> " + exonerate_log)
                    exonerate_fasta_output = "exonerate_results/" + result_subject_name + "/" + query_name + ".fasta"
                    os.system("(tail -n +4 " + exonerate_output + " | head -n -1) > " + exonerate_fasta_output)
