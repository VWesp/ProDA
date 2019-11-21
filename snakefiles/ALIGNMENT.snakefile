import os
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment

rule alignment:
    input:
        "data/queries/{query}.faa",
        "cd_hit/{subject}/{query}_merged.faa"
    output:
        "scores/{subject}/{query}.sc"
    log:
        "log/scores/{subject}/{query}.log"
    run:
        try:
            queries = SeqIO.parse(open(input[0]), "fasta")
            scores = []
            for query in queries:
                scores.append("#" + query.id)
                fasta_sequences = SeqIO.parse(open(input[1]), "fasta")
                for fasta in fasta_sequences:
                    fasta_name = fasta.id.split("::query=")[-1]
                    if(query.id == fasta_name):
                        score = pairwise2.align.globalxx(str(query.seq), str(fasta.seq), score_only=True)
                        scores.append(fasta.id + "\t" + str(score))

                scores.append("\n"*2)

            with open(output[0], "w") as score_writer:
                score_writer.write("\n".join(scores))

            del scores[:]
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(ex))
