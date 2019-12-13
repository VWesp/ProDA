from Bio import Seq

rule Retrieve_Sp_Sequences:
    input:
        "data/subjects/{subject}.fna",
        "alignment/spaln/{subject}/{query}.sp"
    output:
        "alignment/spaln/{subject}/{query}.seq_gff"
    log:
        "log/alignment/spaln/{subject}/{query}.seq_log"
    run:
        try:
            if(os.stat(input[1]).st_size != 0):
                with open(input[1], "r") as sp_reader:
                    subjects = SeqIO.index(input[0], "fasta")
                    sp_results = sp_reader.readlines()
                    contig = None
                    block_start = None
                    query = None
                    orientation = None
                    gff_list = []
                    cds_sequence = []
                    cds_pos = []
                    finished = False
                    for line in sp_results:
                        if(line and not line.startswith(";M")):
                            if(line.startswith(">")):
                                finished = False
                                contig = line.split(" ")[1].split("_cstart:")[0].strip()
                                block_start = int(line.split(" ")[1].split("_cstart:")[1].split("_cend:")[0].strip())
                                query = line.split(" ")[1].split("_query:")[1].strip()
                                orientation = line.split(" ")[2].strip()
                                if(orientation == "-"):
                                    end = int(line.split(" ")[6].strip()) + block_start - 1
                                    start = int(line.split(" ")[8].strip()) + block_start
                                    gff_list.append("#" + contig + "\t" + query + "\t" + str(start) + "\t" + str(end-1) + "\t" + orientation)
                                    gff_list.append(">nuc:" + str(subjects[contig].seq[start:end].reverse_complement()).upper())
                                else:
                                    start = int(line.split(" ")[6].strip()) + block_start - 1
                                    end = int(line.split(" ")[8].strip()) + block_start
                                    gff_list.append("#" + contig + "\t" + query + "\t" + str(start) + "\t" + str(end-1) + "\t" + orientation)
                                    gff_list.append(">nuc:" + str(subjects[contig].seq[start:end]).upper())
                            elif(line.startswith(";C")):
                                if("join" in line):
                                    cds_pos.append(list(filter(None, line.strip().split("join(")[1].split(","))))
                                else:
                                    cds_pos.append(list(filter(None, line.strip().split(";C ")[1].replace(")", "").split(","))))
                            elif(line.startswith("//")):
                                for position in cds_pos:
                                    for start_end_pos in position:
                                        start = int(start_end_pos.split("..")[0].strip()) + block_start - 1
                                        end = int(start_end_pos.split("..")[1].strip()) + block_start
                                        gff_list.append(contig + "\t" + query + "\t" + str(start) + "\t" + str(end) + "\t" + orientation)

                                del cds_pos[:]

                                finished = True
                                gff_list.append(">cds:" + "".join(cds_sequence))
                                gff_list.append(">pep:" + str(Seq("".join(cds_sequence)).translate()))
                                del cds_sequence[:]
                            elif(not finished):
                                cds_sequence.append(line.strip())

                    with open(output[0], "w") as gff_writer:
                        gff_writer.write("\n".join(gff_list))

                    del gff_list[:]

            if(not os.path.exists(output[0])):
                os.system("touch " + output[0])
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()) + "\n" + str(ex))
