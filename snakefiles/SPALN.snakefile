import os
from Bio import SeqIO

rule spaln:
    input:
        "data/queries/{query}.faa",
        "matches/{subject}/{query}.fna"
    output:
        "spaln/{subject}/{query}.sp",
        "spaln/{subject}/{query}.fna",
        temp("spaln/{subject}/{query}_temp.fna"),
        "spaln/{subject}/{query}.faa"
    log:
        "log/spaln/{subject}/{query}.log"
    run:
        fasta_sequences = SeqIO.parse(open(input[1]), "fasta")
        result_list = []
        for fasta in fasta_sequences:
            with open(output[3], "w") as fasta_writer:
                fasta_writer.write(">" + fasta.id + "\n" + str(fasta.seq))

            os.system("(spaln -M -Q3 -O6 -S3 -o" + output[0] + " " + output[3] + " " + input[0] + ") 2>" + log[0])
            with open(output[0], "r") as spaln_reader:
                result_list.append(spaln_reader.read())

        with open(output[1], "w") as spaln_writer:
            spaln_writer.write("\n".join(result_list))

        del result_list[:]
        spaln_results = []
        with open(output[1], "r") as spaln_reader:
            content = spaln_reader.readlines()
            sequence = []
            current_header = None
            current_query = None
            for line in content:
                if(line.startswith(">")):
                    if(current_header != None):
                        spaln_results.append(">" + current_header + "\n" + "".join(sequence))

                    current_header = line.split(" ")[1] + "::query=" + line.split(" ")[10]
                    del sequence[:]
                elif(line.strip().isalpha()):
                    sequence.append(line.strip())

            spaln_results.append(">" + current_header + "\n" + "".join(sequence))
            del sequence[:]

        with open(output[2], "w") as fasta_writer:
            fasta_writer.write("\n".join(spaln_results))

        del spaln_results[:]
        translated_fastas = []
        fastas = SeqIO.parse(open(output[2]), "fasta")
        for fasta in fastas:
            translated_fastas.append(">" + fasta.id + "\n" + str(fasta.seq.translate()))

        with open(output[3], "w") as translated_writer:
            translated_writer.write("\n".join(translated_fastas))

        del translated_fastas[:]
