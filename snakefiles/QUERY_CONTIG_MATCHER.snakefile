import os
from Bio import SeqIO
import traceback

rule Filter_Out_Useless_Regions:
    input:
        "data/subjects/{subject}.fna",
        "blast_results/{subject}/{query}.hit"
    output:
        "matches/{subject}/{query}.fna"
    log:
        "log/matches/{subject}/{query}.log"
    params:
        left=config["left_addendum"],
        right=config["right_addendum"]
    run:
        try:
            if(os.stat(input[1]).st_size != 0):
                subjects_contigs = {}
                subject_name = input[0].split("/")[-1].split(".fna")[0]
                subjects_contigs[subject_name] = {}
                fasta_sequences = SeqIO.parse(open(input[0]), "fasta")
                for fasta in fasta_sequences:
                    subjects_contigs[subject_name][fasta.id] = str(fasta.seq)

                hit_name = input[1].split("/")[-1].split(".hit")[0]
                fasta = []
                with open(input[1], "r") as blast_reader:
                    content = blast_reader.readlines()
                    for line in content:
                        line_splitted = line.split("\t")
                        query_id = line_splitted[0]
                        contig_id = line_splitted[1]
                        hit_start = int(line_splitted[2])
                        hit_end = int(line_splitted[3])

                        if(hit_start > hit_end):
                            temp_loc = hit_start
                            hit_start = hit_end
                            hit_end = temp_loc

                        if(hit_start - params[0] < 0):
                            hit_start = 0
                        else:
                            hit_start -= params[0]

                        if(hit_end + params[1] >= len(subjects_contigs[subject_name][contig_id])):
                            hit_end = len(subjects_contigs[subject_name][contig_id])-1
                        else:
                            hit_end += params[1]

                        header = contig_id + "_start:" + str(hit_start) + "_end:" + str(hit_end) + "_query:" + query_id
                        fasta.append(">" + header + "\n" + sequence)

                subjects_contigs.clear()
                with open(output[0], "w") as match_writer:
                    match_writer.write("\n".join(fasta))

                del fasta[:]

            if(not os.path.exists(output[0])):
                os.system("touch " + output[0])
        except:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()))
