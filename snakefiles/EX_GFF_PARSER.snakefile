from Bio import Seq

rule Retrieve_Ex_Sequences:
    input:
        "data/subjects/{subject}.fna",
        "alignment/exonerate/{subject}/{query}.ryo"
    output:
        "alignment/exonerate/{subject}/{query}.seq_gff"
    log:
        "log/alignment/exonerate/{subject}/{query}.seq_log"
    run:
        try:
            if(os.stat(input[1]).st_size != 0):
                with open(input[1], "r") as ryo_reader:
                    subjects = SeqIO.index(input[0], "fasta")
                    ex_results = ryo_reader.readlines()
                    gff_list = []
                    cds_sequence = []
                    seq_found = False
                    for line in ex_results:
                        if(line and not line.startswith("#")):
                            if(line.startswith(">")):
                                seq_found = True
                            elif(seq_found):
                                cds_sequence.append(line.strip())
                            else:
                                splitted_line = line.split("\t")
                                contig = splitted_line[0].split("_cstart:")[0].strip()
                                block_start = int(splitted_line[0].split("_cstart:")[1].split("_cend:")[0].strip())
                                query = splitted_line[0].split("_query:")[1].strip()
                                type = splitted_line[2]
                                start = int(splitted_line[3].strip()) + block_start - 1
                                end = int(splitted_line[4].strip()) + block_start
                                orientation = splitted_line[6].strip()
                                if(type == "gene"):
                                    gff_list.append("#" + contig + "\t" + query + "\t" + str(start) + "\t" + str(end-1) + "\t" + orientation)
                                    if(orientation == "-"):
                                        gff_list.append(">nuc:" + str(subjects[contig].seq[start:end].reverse_complement()).upper())
                                    else:
                                        gff_list.append(">nuc:" + str(subjects[contig].seq[start:end]).upper())
                                elif(type == "cds"):
                                    gff_list.append(contig + "\t" + query + "\t" + str(start) + "\t" + str(end-1) + "\t" + orientation)
                        elif(line.startswith("#")):
                            seq_found = False
                            if(len(cds_sequence)):
                                gff_list.append(">cds:" + "".join(cds_sequence))
                                gff_list.append(">pep:" + str(Seq("".join(cds_sequence)).translate()))
                                del cds_sequence[:]

                    if(len(cds_sequence)):
                        gff_list.append(">cds:" + "".join(cds_sequence))
                        gff_list.append(">pep:" + str(Seq("".join(cds_sequence)).translate()))
                        del cds_sequence[:]

                    with open(output[0], "w") as gff_writer:
                        gff_writer.write("\n".join(gff_list))

                    del gff_list[:]

            if(not os.path.exists(output[0])):
                os.system("touch " + output[0])
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()) + "\n" + str(ex))
