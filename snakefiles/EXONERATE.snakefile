import os
from Bio import SeqIO

rule exonerate:
    input:
        "data/queries/{query}.faa",
        "matches/{subject}/{query}.fna"
    output:
        "exonerate/{subject}/{query}.ryo",
        "exonerate/{subject}/{query}.faa",
        temp("exonerate/{subject}/{query}.fna")
    params:
        per=config["exonerate_percentage"]
    log:
        "log/exonerate/{subject}/{query}.log"
    run:
        try:
            fastas = []
            if(os.stat(input[1]).st_size != 0):
                #ryo: ::orientation=%g::score=%s::similarity=%ps
                os.system("(exonerate --model protein2genome --targettype dna --querytype protein "
                          "--ryo '>%ti::query=%qi\n%tcs' --showalignment no --showvulgar no "
                          "--refine region --percent " + str(params[0]) + " "
                          "--query " + input[0] + " --target " + input[1] +
                          " > " + output[0] + ") 2> " + log[0])
                os.system("(tail -n +4 " + output[0] + " | head -n -1) > " + output[2])
                fastas = SeqIO.parse(open(output[2]), "fasta")

            translated_fastas = []
            for fasta in fastas:
                translated_fastas.append(">" + fasta.id + "\n" + str(fasta.seq.translate()))

            with open(output[1], "w") as translated_writer:
                translated_writer.write("\n".join(translated_fastas))

            del translated_fastas[:]
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(ex))
