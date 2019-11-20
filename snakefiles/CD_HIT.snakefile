rule cd_hit_exonerate:
    input:
        expand("exonerate/{subject}/{query}.faa", subject=config["subjects"], query=config["queries"])
    output:
        "cd_hit/exonerate/{subject}/{query}.faa"
    log:
        "log/cd_hit/exonerate/{subject}/{query}.log"
    params:
        id=config["cd_hit_threshold"]
    run:
        for fasta_file in input:
            subject_name = fasta_file.split("/")[-2]
            query_name = fasta_file.split("/")[-1].split(".faa")[0]
            cd_hit_output = "cd_hit/exonerate/" + subject_name + "/" + query_name + ".faa"
            cd_hit_log = "log/cd_hit/exonerate/" + subject_name + "/" + query_name + ".log"
            os.system("(cd-hit -i " + fasta_file + " -o " + cd_hit_output + " -c " +
                      str(float(params[0])/100) + ") 2> " + log[0])

rule cd_hit_spaln:
    input:
        expand("spaln/{subject}/{query}.faa", subject=config["subjects"], query=config["queries"])
    output:
        "cd_hit/spaln/{subject}/{query}.faa"
    log:
        "log/cd_hit/spaln/{subject}/{query}.log"
    params:
        id=config["cd_hit_threshold"]
    run:
        for fasta_file in input:
            subject_name = fasta_file.split("/")[-2]
            query_name = fasta_file.split("/")[-1].split(".faa")[0]
            cd_hit_output = "cd_hit/spaln/" + subject_name + "/" + query_name + ".faa"
            cd_hit_log = "log/cd_hit/spaln/" + subject_name + "/" + query_name + ".log"
            os.system("(cd-hit -i " + fasta_file + " -o " + cd_hit_output + " -c " +
                      str(float(params[0])/100) + ") 2> " + log[0])
