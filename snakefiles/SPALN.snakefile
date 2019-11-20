import os
from Bio import SeqIO

rule spaln:
    input:
        queries = expand("data/queries/{query}.faa", query=config["queries"]),
        blast_results = expand("matches/{subject}/{query}.fna", subject=config["subjects"], query=config["queries"])
    output:
        "spaln/{subject}/{query}.faa"
    log:
        "log/spaln/{subject}/{query}.log"
    run:
        result_list = []
        for result in input.blast_results:
            result_subject_name = result.split("/")[-2]
            result_query_name = result.split("/")[-1].split(".fna")[0]
            fasta_sequences = SeqIO.parse(open(result), "fasta")
            for query in input.queries:
                del result_list[:]
                query_name = query.split("/")[-1].split(".faa")[0]
                if(result_query_name == query_name):
                    spaln_output = "spaln/" + "/" + result_subject_name + "/" + query_name + ".sp"
                    spaln_fasta_output = "spaln/" + "/" + result_subject_name + "/" + query_name + ".fna"
                    temp_output = "spaln/" + "/" + result_subject_name + "/" + query_name + "_temp.fna"
                    if(not os.path.exists("spaln/" + result_subject_name)):
                        os.makedirs("spaln/" + result_subject_name)

                    if(not os.path.exists("log/spaln/" + result_subject_name)):
                        os.makedirs("log/spaln/" + result_subject_name)
                        
                    for fasta in fasta_sequences:
                        with open(temp_output, "w") as fasta_writer:
                            fasta_writer.write(">" + fasta.id + "\n" + str(fasta.seq))

                        os.system("(spaln -M -Q3 -O6 -S3 -o" + spaln_output + " " + temp_output + " " + query + ") 2>" + log[0])
                        with open(spaln_output, "r") as spaln_reader:
                            result_list.append(spaln_reader.read())

                        os.remove(temp_output)

                    with open(spaln_output, "w") as spaln_writer:
                        spaln_writer.write("\n".join(result_list))

                    spaln_results = []
                    with open(spaln_output, "r") as spaln_reader:
                        content = spaln_reader.readlines()
                        sequence = []
                        current_header = None
                        for line in content:
                            if(line.startswith(">")):
                                if(current_header != None):
                                    spaln_results.append(">" + current_header + "\n" + "".join(sequence))

                                current_header = line.split(" ")[1] + "::query=" + result_query_name
                                del sequence[:]
                            elif(line.strip().isalpha()):
                                sequence.append(line.strip())

                        spaln_results.append(">" + current_header + "\n" + "".join(sequence))
                        del sequence[:]

                    with open(spaln_fasta_output, "w") as fasta_writer:
                        fasta_writer.write("\n".join(spaln_results))

                    del spaln_results[:]
                    translated_fastas = []
                    fastas = SeqIO.parse(open(spaln_fasta_output), "fasta")
                    for fasta in fastas:
                        translated_fastas.append(">" + fasta.id + "\n" + str(fasta.seq.translate()))

                    translated_spaln_fasta_output = "spaln/" + result_subject_name + "/" + query_name + ".faa"
                    with open(translated_spaln_fasta_output, "w") as translated_writer:
                        translated_writer.write("\n".join(translated_fastas))

                    del translated_fastas[:]
