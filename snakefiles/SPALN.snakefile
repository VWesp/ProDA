import os
from Bio.Seq import Seq

rule spaln:
    input:
        "data/queries/{query}.faa",
        "matches/{subject}/{query}.fna"
    output:
        "spaln/{subject}/{query}.sp",
        "spaln/{subject}/{query}.faa",
        temp("spaln/{subject}/{query}_temp.sp"),
        temp("spaln/{subject}/{query}.fna")
    params:
        pam=config["pam"]
    log:
        "log/spaln/{subject}/{query}.log"
    run:
        try:
            if(os.stat(input[1]).st_size != 0):
                result_list = []
                fasta_sequences = SeqIO.parse(open(input[1]), "fasta")
                for fasta in fasta_sequences:
                    with open(output[3], "w") as fasta_writer:
                        fasta_writer.write(">" + fasta.id + "\n" + str(fasta.seq))

                    os.system("(spaln -M -Q3 -O6 -S3 -yp" + str(params[0]) + " -yq" + str(params[0]) +
                              " -o" + output[2] + " " + output[3] + " " + input[0] + ") 2>" + log[0])
                    with open(output[2], "r") as spaln_reader:
                        result_list.append(spaln_reader.read())

                with open(output[0], "w") as spaln_writer:
                    spaln_writer.write("\n".join(result_list))

                del result_list[:]
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
            else:
                with open(output[0], "w") as empty_writer:
                    empty_writer.write("")

                with open(output[1], "w") as empty_writer:
                    empty_writer.write("")

                with open(output[2], "w") as empty_writer:
                    empty_writer.write("")

                with open(output[3], "w") as empty_writer:
                    empty_writer.write("")
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(ex))
