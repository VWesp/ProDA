rule Filter_Sequences:
    input:
        "merged/{subject}/{query}.seq_gff"
    output:
        "results/above_threshold/{subject}/{query}.gff",
        "results/above_threshold/{subject}/{query}.nuc",
        "results/above_threshold/{subject}/{query}.cds",
        "results/above_threshold/{subject}/{query}.pep",
        "results/below_threshold/{subject}/{query}.gff",
        "results/below_threshold/{subject}/{query}.nuc",
        "results/below_threshold/{subject}/{query}.cds",
        "results/below_threshold/{subject}/{query}.pep"
    log:
        "log/results/above_threshold/{subject}/{query}.log",
        "log/results/below_threshold/{subject}/{query}.log"
    params:
        config["id_threshold"],
        config["sim_threshold"]
    run:
        try:
            if(os.stat(input[0]).st_size != 0):
                above_gff_list = []
                above_nuc_list = []
                above_cds_list = []
                above_pep_list = []
                below_gff_list = []
                below_nuc_list = []
                below_cds_list = []
                below_pep_list = []
                above_index = 0
                below_index = 0
                with open(input[0], "r") as gff_reader:
                    content = gff_reader.readlines()
                    for line in content:
                        if(line):
                            stripped_line = line.strip().split("\t")
                            contig = stripped_line[0]
                            query = stripped_line[1]
                            start = stripped_line[2]
                            end = stripped_line[3]
                            orientation = stripped_line[4]
                            positions = stripped_line[5].split("|")
                            identity = stripped_line[6]
                            similarity = stripped_line[7]
                            nuc_seq = stripped_line[8]
                            cds_seq = stripped_line[9]
                            pep_seq = stripped_line[10]
                            if(float(identity) >= params[0] and float(similarity) >= params[1]):
                                above_gff_list.append("##target\tquery\ttype\tstart\tend\torientation\tidentity\tsimilarity\tindex")
                                above_gff_list.append(contig + "\t" + query + "\t" + "gene" + "\t" + start + "\t" + end + "\t"
                                                      + orientation + "\t" + identity + "\t" + similarity + "\t" + str(above_index))
                                for pos in positions:
                                    cds_start = pos.split(";")[0]
                                    cds_end = pos.split(";")[1]
                                    above_gff_list.append(contig + "\t" + query + "\t" + "cds" + "\t" + cds_start + "\t" + cds_end + "\t"
                                                          + orientation)

                                above_nuc_list.append(">" + contig + "|" + query + "|" + str(above_index)+ "\n" + nuc_seq)
                                above_cds_list.append(">" + contig + "|" + query + "|" + str(above_index) + "\n" + cds_seq)
                                above_pep_list.append(">" + contig + "|" + query + "|" + str(above_index) + "\n" + pep_seq)
                                above_index += 1
                            else:
                                below_gff_list.append("##target\tquery\ttype\tstart\tend\torientation\tidentity\tsimilarity\tindex")
                                below_gff_list.append(contig + "\t" + query + "\t" + "gene" + "\t" + start + "\t" + end + "\t"
                                                      + orientation + "\t" + identity + "\t" + similarity + "\t" + str(below_index))
                                for pos in positions:
                                    cds_start = pos.split(";")[0]
                                    cds_end = pos.split(";")[1]
                                    below_gff_list.append(contig + "\t" + query + "\t" + "cds" + "\t" + cds_start + "\t" + cds_end + "\t"
                                                          + orientation)

                                below_nuc_list.append(">" + contig + "|" + query + "|" + str(below_index) + "\n" + nuc_seq)
                                below_cds_list.append(">" + contig + "|" + query + "|" + str(below_index) + "\n" + cds_seq)
                                below_pep_list.append(">" + contig + "|" + query + "|" + str(below_index) + "\n" + pep_seq)
                                below_index += 1

                with open(output[0], "w") as gff_writer:
                    gff_writer.write("\n".join(above_gff_list))

                with open(output[1], "w") as nuc_writer:
                    nuc_writer.write("\n".join(above_nuc_list))

                with open(output[2], "w") as cds_writer:
                    cds_writer.write("\n".join(above_cds_list))

                with open(output[3], "w") as pep_writer:
                    pep_writer.write("\n".join(above_pep_list))

                with open(output[4], "w") as gff_writer:
                    gff_writer.write("\n".join(below_gff_list))

                with open(output[5], "w") as nuc_writer:
                    nuc_writer.write("\n".join(below_nuc_list))

                with open(output[6], "w") as cds_writer:
                    cds_writer.write("\n".join(below_cds_list))

                with open(output[7], "w") as pep_writer:
                    pep_writer.write("\n".join(below_pep_list))

                del above_gff_list[:]
                del above_nuc_list[:]
                del above_cds_list[:]
                del above_pep_list[:]
                del below_nuc_list[:]
                del below_cds_list[:]
                del below_pep_list[:]

            for out in output:
                if(not os.path.exists(out)):
                    os.system("touch " + out)
        except Exception as ex:
            with open(log[0], "w") as log_writer:
                log_writer.write(str(traceback.format_exc()) + "\n" + str(ex))
