import os

rule cd_hit:
    input:
        "exonerate/{subject}/{query}.faa",
        "spaln/{subject}/{query}.faa"
    output:
        "cd_hit/{subject}/{query}_joined.faa",
        "cd_hit/{subject}/{query}_merged.faa",
        temp("cd_hit/{subject}/{query}_merged.clstr")
    log:
        "log/cd_hit/{subject}/{query}.log"
    params:
        id=config["cd_hit_threshold"]
    run:
        try:
            merge = []
            with open(input[0], "r") as input_reader:
                merge.append(input_reader.read())

            with open(input[1], "r") as input_reader:
                merge.append(input_reader.read())

            with open(output[0], "w") as input_writer:
                input_writer.write("\n".join(merge))

            del merge[:]
            os.system("(cd-hit -i " + output[0] + " -o " + output[1] +
                      " -c " + str(float(params[0])/100) + ") 2> " + log[0])
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(ex))
