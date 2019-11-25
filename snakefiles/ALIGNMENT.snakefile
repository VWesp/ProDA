import os
from Bio import SeqIO

rule alignment:
    input:
        "data/queries/{query}.faa",
        "cd_hit/{subject}/{query}_merged.faa"
    output:
        "scores/{subject}/{query}.sim",
        temp("scores/{subject}/{query}_temp.st"),
        temp("scores/{subject}/{query}_temp_query.faa"),
        temp("scores/{subject}/{query}_temp_hit.faa")
    log:
        "log/scores/{subject}/{query}.log"
    run:
        try:
            if(os.stat(input[1]).st_size != 0):
                scores = []
                queries = SeqIO.parse(open(input[0]), "fasta")
                for query in queries:
                    scores.append("#" + query.id)
                    fasta_sequences = SeqIO.parse(open(input[1]), "fasta")
                    for fasta in fasta_sequences:
                        fasta_name = fasta.id.split("::query=")[-1]
                        if(query.id == fasta_name):
                            with open(output[2], "w") as query_writer:
                                query_writer.write(">" + query.id + "\n" + str(query.seq))

                            with open(output[3], "w") as hit_writer:
                                hit_writer.write(">" + fasta.id + "\n" + str(fasta.seq))

                            os.system("(stretcher -asequence " + output[2] + " -sprotein1 -bsequence " +
                                      output[3] + " -auto -stdout > " + output[1] + ") 2> " + log[0])
                            with open(output[1], "r") as sim_reader:
                                content = sim_reader.readlines()
                                for line in content:
                                    if(line.startswith("# Similarity")):
                                        similarity = line.split("(")[1][:-3]
                                        scores.append(fasta.id + "\t" + similarity + "\t" + str(fasta.seq) + "\t" + str(query.seq))
                                        break

                    scores.append("\n"*2)

                with open(output[0], "w") as score_writer:
                    score_writer.write("\n".join(scores))

                del scores[:]
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
